using CryspBrainSim
using Test
using NPZ
using LinearAlgebra: diag
using SpecialFunctions: erfc

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
end
