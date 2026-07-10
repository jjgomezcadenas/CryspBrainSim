using CryspBrainSim
using Test
using NPZ
using HDF5
using TOML
using LinearAlgebra: diag
using SpecialFunctions: erfc, erfcinv

# Cross-validation reference (dev/PLAN.md validation ladder rung 3): shared
# numpy-generated arrays plus the frozen Python fit outputs, produced once by
# test/make_reference.py from the frozen py/ estimator and committed.
const REF = npzread(joinpath(@__DIR__, "data", "endpoint_reference.npz"))

# Documented tolerances, measured on 2026-07-05 (Julia 1.12 / LsqFit 0.16.1 vs
# scipy 1.17.1) and set an order of magnitude above the measured agreement:
#   depth_profile         bit-exact on integer-count images   → atol 1e-9
#   fit params / R levels measured ≤ 2e-5 relative            → rtol 1e-3
#   z0_err                measured ≤ 5e-5 relative            → rtol 1e-3
#   pcov diagonal         measured ≤ 1e-4 relative            → rtol 1e-2
#   sigma_R               exact arithmetic                     → rtol 1e-12
# The two LM minimisers stop at slightly different points; the ~1e-5 parameter
# scatter is that stopping difference, not a convention mismatch.
const RTOL_FIT = 1e-3
const RTOL_PCOV = 1e-2

erfc_truth(z; base, amp, z0, w) = @. base + amp * 0.5 * erfc((z - z0) / (sqrt(2.0) * w))

# ---------------------------------------------------------------------------
# Synthetic products fixtures (step-3 tests run without PtCryspProds on disk)
# ---------------------------------------------------------------------------

# A schema-faithful miniature lors_shardNNN.h5: 8 LORs on a 387 mm ring —
# 5 true (row 1 degenerate: identical endpoints), 1 single scatter,
# 1 multiple scatter, 1 random; one LOR (row 8) outside the ±3 ns window.
function write_mini_shard(path; realization=0, crystal="BGO", nevents=40,
                          scenario="mini_scen", drop_attr=nothing)
    n = 8
    q(v) = Int16.(round.(v ./ 0.1))                      # quantize mm → Int16
    x1 = q([390.0, 390, 390, 390, 390, 390, 390, 390])
    y1 = q(zeros(n))
    z1 = q([10.0, 20, 30, 40, 50, 60, 70, 80])
    x2 = q([390.0, -390, -390, -390, -390, -390, -390, -390])  # row 1 degenerate
    y2 = q(zeros(n))
    z2 = copy(z1); z2[1] = z1[1]                          # row 1: same endpoint
    truth = Int8[0, 0, 0, 0, 1, 1, 2, 0]
    nscat1 = Int8[0, 0, 0, 0, 1, 1, 0, 0]
    nscat2 = Int8[0, 0, 0, 0, 0, 1, 0, 0]
    e = Int16.(round.([511.0, 511, 500, 480, 460, 455, 470, 511] ./ 0.1))
    t1 = Float32[1.0, 1, 1, 1, 1, 1, 1, 6]
    t2 = Float32[0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]  # row 8: Δt = 5.5 ns
    x0 = q([0.0, 1, -1, 2, -2, 3, -3, 4])
    y0 = q([0.0, -5, 5, -10, 10, -15, 15, -20])
    z0 = q([-90.0, -60, -30, 0, 30, 60, 90, 95])
    h5open(path, "w") do f
        for (k, v) in (
            "event" => Int32.(1:n), "truth" => truth,
            "x1_mm" => x1, "y1_mm" => y1, "z1_mm" => z1,
            "e1_keV" => e, "t1_ns" => t1,
            "iz1" => Int16.(ones(n)), "iphi1" => Int16.(ones(n)), "nscat1" => nscat1,
            "x2_mm" => x2, "y2_mm" => y2, "z2_mm" => z2,
            "e2_keV" => reverse(e), "t2_ns" => t2,
            "iz2" => Int16.(ones(n)), "iphi2" => Int16.(2 .* ones(n)), "nscat2" => nscat2,
            "dt_ns" => Float32.(t1 .- t2), "x0_mm" => x0, "y0_mm" => y0, "z0_mm" => z0)
            write(f, k, v)
        end
        a = Dict{String,Any}(
            "scenario" => scenario, "crystal" => crystal, "budget" => "fast",
            "dose_Gy" => 1.0, "master_seed" => 1, "realization" => realization,
            "n_phi" => 48, "n_z" => 20, "tau_ns" => 3.0, "emin_keV" => 450.0,
            "eres" => 0.1, "sigma_xyz_mm" => 1.7, "nevents" => nevents,
            "nrows" => n, "xyz_scale_mm" => 0.1, "e_scale_keV" => 0.1)
        drop_attr !== nothing && delete!(a, drop_attr)
        for (k, v) in a
            HDF5.attributes(f)[k] = v
        end
    end
    return path
end

# A miniature truth/ bundle with analytic edges: dose erfc at 150 mm (w 2),
# activity erfc at 140 mm (w 3) — activity-R50 = 140 exactly by construction.
function write_mini_truth(scenario_dir; z_shift_activity=0.0)
    tdir = joinpath(scenario_dir, "truth")
    mkpath(tdir)
    z = collect(0.0:1.0:200.0)
    dose = @. 0.5 * erfc((z - 150.0) / (sqrt(2.0) * 2.0))
    act = @. 1000.0 * 0.5 * erfc((z - 140.0) / (sqrt(2.0) * 3.0))
    open(joinpath(tdir, "depth_dose.csv"), "w") do io
        println(io, "z_mm,edep_total_MeV,edep_primary_MeV,edep_core_MeV,dose_core_Gy")
        for i in eachindex(z)
            println(io, "$(z[i]),$(1e6 * dose[i]),$(9e5 * dose[i]),$(1e3 * dose[i]),$(dose[i])")
        end
    end
    open(joinpath(tdir, "activity_profile_fast.csv"), "w") do io
        println(io, "z_mm,O15,C11,total")
        for i in eachindex(z)
            za = z[i] + z_shift_activity
            println(io, "$(za),$(0.7 * act[i]),$(0.3 * act[i]),$(act[i])")
        end
    end
    return scenario_dir
end

@testset "CryspBrainSim" begin
    @testset "package" begin
        @test isdefined(CryspBrainSim, :RecoCryspTools)
    end

    @testset "depth_profile — cross-validation on the shared image" begin
        img = Float64.(REF["image"])
        z, prof = depth_profile(img;
                                voxel_size_mm=Tuple(REF["image_vox"]),
                                beam_axis=3,
                                roi_radius_mm=REF["image_roi_radius"][1])
        @test z ≈ REF["zB"] atol = 1e-9
        @test prof ≈ REF["profB"] atol = 1e-9

        # The disc ROI keeps a strict subset of the whole-plane sum.
        _, prof_all = depth_profile(img; voxel_size_mm=Tuple(REF["image_vox"]),
                                    beam_axis=3)
        @test all(prof .<= prof_all)
        @test sum(prof) < sum(prof_all)

        # beam_axis handling: permuting the image and the axis label together
        # leaves the profile unchanged.
        imgp = permutedims(img, (3, 1, 2))
        _, profp = depth_profile(imgp;
                                 voxel_size_mm=(1.5, 1.5, 1.5), beam_axis=1,
                                 roi_radius_mm=REF["image_roi_radius"][1])
        @test profp == prof

        @test_throws ArgumentError depth_profile(img; voxel_size_mm=(1.5, 1.5),
                                                 beam_axis=3)
        @test_throws ArgumentError depth_profile(img;
                                                 voxel_size_mm=(1.5, 1.5, 1.5),
                                                 beam_axis=4)
    end

    @testset "distal_window" begin
        @test distal_window(-5.0) == (-25.0, 10.0)
        @test distal_window(150.0; proximal_margin_mm=25.0,
                            distal_margin_mm=20.0) == (125.0, 170.0)
    end

    @testset "fit_endpoint — 1-D edge (test A), weighted" begin
        win = (REF["winA"][1], REF["winA"][2])
        f = fit_endpoint(REF["zA"], REF["profA"]; window=win, weighted=true)
        @test f.n_points == REF["A_weighted_n_points"][1]
        @test f.popt ≈ REF["A_weighted_popt"] rtol = RTOL_FIT
        @test f.z0_err ≈ REF["A_weighted_z0_err"][1] rtol = RTOL_FIT
        @test diag(f.pcov) ≈ REF["A_weighted_pcov_diag"] rtol = RTOL_PCOV
        @test [f.R[0.5], f.R[0.8], f.R[0.2]] ≈ REF["A_weighted_R"] rtol = RTOL_FIT
        @test f.R[0.5] == f.z0                       # erfcinv(1) = 0
        @test f.R[0.8] < f.z0 < f.R[0.2]             # level ordering

        # Truth recovery: the data were generated with z0 = 150, w = 3.
        @test abs(f.z0 - 150.0) < 3 * f.z0_err
        @test abs(f.w - 3.0) < 0.5
    end

    @testset "fit_endpoint — 1-D edge (test A), unweighted covariance" begin
        # scipy default (absolute_sigma=False, MSE-scaled pcov) vs LsqFit's
        # unweighted vcov — the second half of the covariance-convention check.
        win = (REF["winA"][1], REF["winA"][2])
        f = fit_endpoint(REF["zA"], REF["profA"]; window=win, weighted=false)
        @test f.popt ≈ REF["A_unweighted_popt"] rtol = RTOL_FIT
        @test f.z0_err ≈ REF["A_unweighted_z0_err"][1] rtol = RTOL_FIT
    end

    @testset "fit_endpoint — image profile (test B), weighted" begin
        img = Float64.(REF["image"])
        z, prof = depth_profile(img;
                                voxel_size_mm=Tuple(REF["image_vox"]),
                                beam_axis=3,
                                roi_radius_mm=REF["image_roi_radius"][1])
        f = fit_endpoint(z, prof; window=(REF["winB"][1], REF["winB"][2]),
                         weighted=true)
        @test f.n_points == REF["B_weighted_n_points"][1]
        @test f.popt ≈ REF["B_weighted_popt"] rtol = RTOL_FIT
        @test f.z0_err ≈ REF["B_weighted_z0_err"][1] rtol = RTOL_FIT
        @test [f.R[0.5], f.R[0.8], f.R[0.2]] ≈ REF["B_weighted_R"] rtol = RTOL_FIT
        @test abs(f.z0 - 150.0) < 1.0                # sub-voxel truth recovery
        @test 0.0 < f.z0_err < 0.5                   # sub-voxel fit error
    end

    @testset "fit_endpoint — failure paths" begin
        # Fewer than 4 points inside the window → NaN endpoints, count kept.
        f = fit_endpoint(REF["zA"], REF["profA"]; window=(0.0, 1.0))
        @test isnan(f.z0) && isnan(f.w) && isnan(f.z0_err)
        @test all(isnan, values(f.R))
        @test f.n_points == 2
        @test f.popt === nothing && f.pcov === nothing
    end

    @testset "sigma_R" begin
        s = sigma_R(REF["sig_endpoints"];
                    dose_bragg_peak=REF["sig_dose_bragg_peak"][1])
        @test [s.n_ok, s.n_fail, s.mean, s.sigma, s.sem, s.offset] ≈
              REF["sig_ref"] rtol = 1e-12

        @test sigma_R([150.0, 151.0]).offset === nothing

        few = sigma_R([150.0, NaN, Inf])
        @test few.n_ok == 1 && few.n_fail == 2
        @test isnan(few.mean) && isnan(few.sigma) && isnan(few.sem)
        @test few.offset === nothing
    end

    @testset "read_csv_table" begin
        mktempdir() do dir
            p = joinpath(dir, "t.csv")
            write(p, "name,a,b\nfoo,1.5,2\nbar,-3e2,4\n")
            t = read_csv_table(p)
            @test t["name"] == ["foo", "bar"]
            @test t["a"] == [1.5, -300.0]
            @test t["b"] == [2.0, 4.0]
            write(p, "a,b\n1,2\n3\n")
            @test_throws ErrorException read_csv_table(p)
        end
    end

    @testset "products — shards on a miniature tree" begin
        mktempdir() do dir
            leaf = joinpath(dir, "mini_scen", "ring", "bgo", "fast_1Gy")
            mkpath(leaf)
            s0 = write_mini_shard(joinpath(leaf, "lors_shard000.h5"); realization=0)
            s1 = write_mini_shard(joinpath(leaf, "lors_shard001.h5"); realization=1)

            @test leaf_dir(dir; scenario="mini_scen", scanner="ring",
                           crystal="bgo", leaf="fast_1Gy") == leaf
            @test_throws ErrorException leaf_dir(dir; scenario="mini_scen",
                                                 scanner="ring", crystal="csi",
                                                 leaf="fast_1Gy")
            @test shard_files(leaf) == [s0, s1]

            a = shard_attrs(s0)
            @test a["realization"] == 0 && a["crystal"] == "BGO"
            bad = write_mini_shard(joinpath(dir, "bad.h5"); drop_attr="master_seed")
            @test_throws ErrorException shard_attrs(bad)

            r = read_shard(s0)
            @test r.n_dropped == 1                    # the degenerate row
            @test length(r.coinc) == 7

            pool = pool_shards([s0, s1])
            @test length(pool.coinc) == 14
            @test pool.n_dropped == 2
            @test length(pool.shard_attrs) == 2

            # Not one master: same realization index twice / common-mode drift.
            @test_throws ErrorException pool_shards([s0, s0])
            s2 = write_mini_shard(joinpath(leaf, "lors_shard002.h5");
                                  realization=2, crystal="CsI")
            @test_throws ErrorException pool_shards([s0, s1, s2])
        end
    end

    @testset "products — geometry, phantom, truth_dir" begin
        mktempdir() do dir
            gj = joinpath(dir, "scanner_geometry.json")
            write(gj, """{"scanner": {"name": "ring", "shape": "cyl_shell",
                  "r_inner_cm": 38.7, "wall_thickness_cm": 3.7,
                  "half_length_cm": 51.2, "n_phi": 48, "n_z": 20}}""")
            g = scanner_geometry(dir)
            @test g.name == "ring" && g.shape == "cyl_shell"
            @test g.r_inner_mm == 387.0 && g.wall_mm == 37.0
            @test g.half_length_mm == 512.0 && g.n_phi == 48 && g.n_z == 20

            ph = joinpath(dir, "phantom"); mkpath(ph)
            write(joinpath(ph, "phantom_regions.csv"),
                  "region,priority,material,solid,a_mm,b_mm,c_mm,cx_mm,cy_mm,cz_mm," *
                  "euler_x_deg,euler_y_deg,euler_z_deg\n" *
                  "head,0,G4_BRAIN_ICRP,ellipsoid,72,87,102,0,-30,0,0,0,0\n")
            write(joinpath(ph, "phantom_material_g4_brain_icrp_meta.csv"),
                  "material,energy_keV,density_g_cm3,mean_excitation_eV,mu_rho_cm2_g," *
                  "mu_cm_inv,mu_mm_inv,mean_free_path_cm,note\n" *
                  "G4_BRAIN_ICRP,511.0,1.04,73.3,0.0953,0.0991,0.009913,10.09,note\n")
            reg = phantom_region(dir)
            @test reg.material == "G4_BRAIN_ICRP" && reg.solid == "ellipsoid"
            @test reg.semi_axes == (72.0, 87.0, 102.0) && reg.centre == (0.0, -30.0, 0.0)
            mu = material_mu(dir, "G4_BRAIN_ICRP")
            @test mu.mu_mm_inv ≈ 0.009913 && mu.energy_keV == 511.0

            @test_throws ErrorException truth_dir(dir)
            mkpath(joinpath(dir, "truth"))
            @test truth_dir(dir) == joinpath(dir, "truth")
        end
    end

    @testset "shard_stats — ShardStats on the miniature shard" begin
        mktempdir() do dir
            f = write_mini_shard(joinpath(dir, "lors_shard000.h5"); nevents=40)
            q = shard_stats(f; r_inner_mm=387.0)
            @test q.nrows == 8 && q.nevents == 40
            @test q.acceptance ≈ 0.2
            @test q.n_true == 5 && q.n_random == 1
            @test q.n_scatter_single == 1 && q.n_scatter_multiple == 1
            @test q.frac_true ≈ 5 / 8 && q.frac_scatter ≈ 2 / 8 && q.frac_random ≈ 1 / 8
            @test q.n_degenerate == 1
            @test q.e_range == (455.0, 511.0)
            @test q.r_range == (390.0, 390.0)
            @test q.doi_range == (3.0, 3.0)
            @test q.dt_frac_in_tau ≈ 7 / 8       # row 8 sits outside ±3 ns
            @test q.tau_ns == 3.0
            @test q.src_min == (-3.0, -20.0, -90.0)
            @test q.src_max == (4.0, 15.0, 95.0)
        end
    end

    @testset "characterize — analytic truth bundle" begin
        mktempdir() do dir
            write_mini_truth(dir)

            # Analytic expectations on the 1 mm grid: dose-R80 solves
            # 0.5·erfc(u) = 0.8 at edge (150, w=2); activity-R50 is exactly
            # the edge position 140 (erfc(0)/2 = 1/2). Linear interpolation on
            # the grid reads them to a few 0.01 mm.
            R80_exact = 150.0 + sqrt(2.0) * 2.0 * erfcinv(1.6)
            ref = characterize(dir)
            @test abs(ref.dose_R80 - R80_exact) < 0.05
            @test abs(ref.activity_R50 - 140.0) < 0.05
            @test abs(ref.activity_R50_fit - 140.0) < 0.05
            @test abs(ref.activity_w_fit - 3.0) < 0.05
            @test ref.offset ≈ ref.activity_R50 - ref.dose_R80
            @test ref.window == (ref.activity_R50 - 20.0, ref.activity_R50 + 15.0)

            path = write_reference(ref, joinpath(dir, "out", "truth_reference.toml"))
            back = TOML.parsefile(path)
            @test back["dose_R80_mm"] ≈ ref.dose_R80
            @test back["activity_R50_mm"] ≈ ref.activity_R50
            @test back["offset_mm"] ≈ ref.offset
            @test back["window_mm"] ≈ collect(ref.window)
        end

        # Broken bundle: activity on a shifted z-frame errors.
        mktempdir() do dir
            write_mini_truth(dir; z_shift_activity=0.5)
            @test_throws ErrorException characterize(dir)
        end

        # distal_crossing reads the LAST falling edge and NaNs when absent.
        z = collect(0.0:1.0:10.0)
        y = [0.0, 1, 0.2, 1, 1, 1, 1, 1, 0.5, 0.0, 0.0]  # dip + final edge
        @test distal_crossing(z, y; level=0.5) ≈ 8.0
        @test isnan(distal_crossing(z, fill(1.0, 11); level=0.5))
    end

    @testset "mumap — chords, factors, voxel route" begin
        @test centered_grid((41, 41, 160), (1.5, 1.5, 1.5)) ==
              (-30.0f0, -30.0f0, -119.25f0)

        axes3 = (72.0f0, 87.0f0, 102.0f0)
        ctr = (0.0f0, -30.0f0, 0.0f0)
        μ = 0.009913f0

        # Axial LOR through the ellipsoid centre: chord = 2c = 204 mm.
        xs = reshape(Float32[0, -30, -500], 3, 1)
        xe = reshape(Float32[0, -30, 500], 3, 1)
        a = attenuation_ellipsoid(xs, xe; semi_axes=axes3, centre=ctr, mu_mm_inv=μ)
        @test a[1] ≈ exp(-μ * 204.0f0) rtol = 1e-5

        # A LOR clear of the body survives untouched.
        miss = attenuation_ellipsoid(reshape(Float32[300, 300, -500], 3, 1),
                                     reshape(Float32[300, 300, 500], 3, 1);
                                     semi_axes=axes3, centre=ctr, mu_mm_inv=μ)
        @test miss[1] == 1.0f0

        # Voxel route agrees with the analytic route to voxelization accuracy.
        n = (80, 120, 105)
        vs = (2.0f0, 2.0f0, 2.0f0)
        org = centered_grid(n, vs)
        mumap = build_mumap(; n=n, img_origin=org, voxsize=vs,
                            semi_axes=axes3, centre=ctr, mu_mm_inv=μ)
        vol_vox = sum(mumap .> 0) * prod(vs)
        vol_true = 4π / 3 * prod(axes3)
        @test abs(vol_vox / vol_true - 1) < 0.01

        many_s = Float32[0 -30 -500; 20 0 -500; 0 -100 -500; 60 -30 -400]'
        many_e = Float32[0 -30 500; -20 -60 500; 0 40 500; -60 -30 400]'
        av = attenuation_mumap(Float32.(many_s), Float32.(many_e), mumap, org, vs)
        ae = attenuation_ellipsoid(Float32.(many_s), Float32.(many_e);
                                   semi_axes=axes3, centre=ctr, mu_mm_inv=μ)
        @test all(abs.(log.(av) .- log.(ae)) .< μ * 2 * maximum(vs))  # ≤ 2-voxel chord error
    end

    @testset "config — the frozen run parameters load typed" begin
        k = load_run_parameters()
        @test k.grid.n == (64, 64, 96)
        @test k.grid.img_origin == (-47.25f0, -47.25f0, -119.25f0)
        @test k.grid.voxsize == (1.5f0, 1.5f0, 1.5f0)
        # No radius: the profile read is whole-plane (the settled protocol).
        @test k.roi.radius_mm === nothing && k.roi.centre_mm == (0.0, 0.0)
        @test k.window[1] < k.window[2] < 0
        @test k.niter == 50
        @test k.n_sens == 1_000_000_000
        @test k.truth_selection == "trues-only"
    end

    @testset "thinning — Bernoulli mask, seed namespace, dose anchor" begin
        # Bit-for-bit reproducibility per realization index; independence
        # across indices.
        m1 = thin_mask(100_000, 10_000, 7)
        @test m1 == thin_mask(100_000, 10_000, 7)
        m2 = thin_mask(100_000, 10_000, 8)
        @test m1 != m2

        # Binomial count statistics: kept ~ Binomial(M, p); 5σ tolerance.
        M, target = 400_000, 40_000
        σ = sqrt(target * (1 - target / M))
        counts = [count(thin_mask(M, target, k)) for k in 1:20]
        @test all(abs.(counts .- target) .< 5σ)
        # The count FLUCTUATES across realizations (Bernoulli, not exact-count).
        @test length(unique(counts)) > 1
        # The ensemble mean sits on the target (sem = σ/√20).
        @test abs(sum(counts) / 20 - target) < 5σ / sqrt(20)

        # Two realizations overlap at rate ≈ p² (independent draws).
        both = count(m1 .& m2)
        @test abs(both - 100_000 * 0.1^2) < 5 * sqrt(100_000 * 0.01)

        # Bounds and errors.
        @test count(thin_mask(1000, 0, 1)) == 0
        @test count(thin_mask(1000, 1000, 1)) == 1000
        @test_throws ArgumentError thin_mask(1000, 1001, 1)

        # The dose anchor: one acquisition's count scaled by the dose ratio.
        @test dose_to_counts(1.0, 1.0, 174_296_897, 10) == 17_429_690
        @test dose_to_counts(0.1, 1.0, 174_296_897, 10) == 1_742_969
        @test_throws ArgumentError dose_to_counts(2.0, 1.0, 100, 10)
        @test_throws ArgumentError dose_to_counts(0.0, 1.0, 100, 10)
    end

    @testset "sensitivity_cache_name + dose_tag" begin
        # The scenario and ring are carried by the path; the name holds only
        # the grid, its origin, and n_sens.
        p = (grid=(n=(64, 64, 96), voxsize=(1.5f0, 1.5f0, 1.5f0),
                   img_origin=(-47.25f0, -47.25f0, -119.25f0)),
             n_sens=1_000_000_000)
        @test sensitivity_cache_name(p) ==
              "grid64x64x96_1.5mm_orgm47.25_m47.25_m119.25_n1000000000"
        # n_sens override for the builder's sweep of the sample size.
        @test endswith(sensitivity_cache_name(p; n_sens=500_000_000), "_n500000000")
        # A different origin yields a different name (no collision).
        p2 = (grid=(n=(64, 64, 96), voxsize=(1.5f0, 1.5f0, 1.5f0),
                    img_origin=(0.0f0, 0.0f0, 0.0f0)), n_sens=1_000_000_000)
        @test sensitivity_cache_name(p2) != sensitivity_cache_name(p)

        @test dose_tag(1.0) == "1Gy"
        @test dose_tag(0.5) == "0p5Gy"
        @test dose_tag(0.05) == "0p05Gy"
        @test dose_tag(2.0) == "2Gy"
    end

    @testset "output layout — crystal_label + path helpers" begin
        @test crystal_label("BGO", 37.0) == "bgo_3X0"   # 37 / 11.18 ≈ 3.3
        @test crystal_label("BGO", 22.4) == "bgo_2X0"   # 22.4 / 11.18 ≈ 2.0
        @test crystal_label("CsI", 37.2) == "csi_2X0"   # 37.2 / 18.6 = 2.0
        @test_throws ErrorException crystal_label("LYSO", 30.0)

        r = "/tmp/out"
        @test scenario_out("scen"; root=r) == joinpath(r, "scen")
        @test truth_out("scen"; root=r) == joinpath(r, "scen", "truth")
        @test sensitivity_out("scen", "closed", "ring1"; root=r) ==
              joinpath(r, "scen", "closed", "ring1", "sensitivity")
        @test config_out("scen", "closed", "ring1", "bgo_3X0"; root=r) ==
              joinpath(r, "scen", "closed", "ring1", "bgo_3X0")
        @test validation_out(; root=r) == joinpath(r, "validation")
    end

    @testset "scanner descriptors — spec + written TOML" begin
        geo = (name="crysp_ring_1m", shape="cyl_shell", r_inner_mm=387.0,
               wall_mm=37.0, half_length_mm=512.0, n_phi=48, n_z=20)
        spec = scanner_spec(geo)
        @test spec["r_outer_mm"] == 424.0
        @test spec["length_mm"] == 1024.0
        @test spec["n_crystals"] == 960
        @test spec["crystal_axial_mm"] ≈ 1024.0 / 20
        @test spec["crystal_transverse_mm"] ≈ 2π * 387.0 / 48
        @test spec["crystal_radial_mm"] == 37.0

        mktempdir() do dir
            gp = write_ring_geometry(geo; scenario="scen", topology="closed",
                                     ring="ring1", root=dir)
            @test gp == joinpath(dir, "scen", "closed", "ring1", "geometry.toml")
            back = TOML.parsefile(gp)
            @test back["n_crystals"] == 960 && back["shape"] == "cyl_shell"

            det = (energy_resolution_fwhm=0.10, sigma_xyz_mm=1.7,
                   emin_keV=450.0, tau_ns=3.0)
            cp = write_crystal_spec(; scenario="scen", topology="closed",
                                    ring="ring1", crystal="bgo_3X0",
                                    material="bgo", wall_mm=37.0, detector=det,
                                    root=dir)
            @test cp == joinpath(dir, "scen", "closed", "ring1", "bgo_3X0",
                                 "crystal.toml")
            cb = TOML.parsefile(cp)
            @test cb["material"] == "BGO" && cb["label"] == "bgo_3X0"
            @test cb["thickness_X0"] ≈ 37.0 / 11.18
            @test cb["detector"]["energy_resolution_fwhm"] == 0.10
            @test cb["detector"]["sigma_xyz_mm"] == 1.7
            @test cb["detector"]["emin_keV"] == 450.0
            @test cb["detector"]["tau_ns"] == 3.0
        end
    end

    @testset "sensitivity — chunked base, scale, cache roundtrip" begin
        n = (16, 16, 16)
        vs = (8.0f0, 8.0f0, 8.0f0)
        org = centered_grid(n, vs)
        atten_ones(xs, xe) = ones(Float32, size(xs, 2))
        args = (r_inner_mm=200.0, half_length_mm=200.0, n=n, img_origin=org,
                voxsize=vs, attenuation=atten_ones, n_sens=200_000,
                chunk=60_000, seed=7, progress=false)
        base = sensitivity_base(; args...)
        @test size(base) == n
        @test all(isfinite, base)
        mid = n .÷ 2
        @test base[mid...] > 0            # the FOV centre is illuminated

        # Same seed reproduces to the atomic-accumulation tolerance
        # (order-nondeterministic Float32 adds — compare with rtol, never ==).
        base2 = sensitivity_base(; args...)
        @test base ≈ base2 rtol = 1e-3

        # The per-realization scale.
        sens = scaled_sensitivity(base, 100_000, 200_000)
        @test sens[mid...] ≈ 0.5f0 * base[mid...]

        mktempdir() do dir
            p = save_sensitivity(joinpath(dir, "sens", "base"), base;
                                 r_inner_mm=200.0, half_length_mm=200.0,
                                 n_sens=200_000, chunk=60_000, seed=7,
                                 img_origin=org, voxsize=vs,
                                 attenuation_meta=(route="none", mu_mm_inv=0.0))
            b2, meta = load_sensitivity(p)
            @test b2 == base
            @test meta["draw"]["n_sens"] == 200_000
            @test meta["grid"]["n"] == collect(n)
            @test meta["scanner"]["r_inner_mm"] == 200.0
            @test meta["attenuation"]["route"] == "none"
            @test haskey(meta, "recocrysp_sha")
        end
    end
end
