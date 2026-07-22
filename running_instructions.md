# Running the statistical-procedure jobs


Use `drivers/statistical_procedure_jobs.jl`. It writes one durable result after
every completed reconstruction. Do not use the older monolithic
`drivers/statistical_procedure.jl` for production runs.

## GPU strategy

Use one GPU MLEM reconstruction at a time:

```bash
julia -t 1 --project=. ...
```

Do not launch several GPU MLEM jobs concurrently. Existing output files are
skipped automatically, so a stopped job can be restarted safely.

## 1. Independent 1-Gy shards

Run the ten shards sequentially:

```bash
for i in {0..9}; do
  julia -t 1 --project=. drivers/statistical_procedure_jobs.jl \
    --stage shard --index "$i"
done
```

For each completed shard the driver writes:

```text
out/.../statistical_procedure/del120s_ac300s_1Gy_D1p0Gy/shards/
    shardNNN.toml
    shardNNN.png
```

The TOML contains event count, fitted `R50`, fit uncertainty, `chi2/dof`, erfc
parameters, profile, and fit window. The PNG shows the profile and fitted erfc
edge.

## 2. Washed thinning ensemble

Run washed realizations in manageable sequential batches. The first ten of an
`N=100` ensemble are:

```bash
julia -t 1 --project=. drivers/statistical_procedure_jobs.jl \
  --stage ensemble --mode washed --first 1 --last 10
```

Continue with successive ranges:

```bash
julia -t 1 --project=. drivers/statistical_procedure_jobs.jl \
  --stage ensemble --mode washed --first 11 --last 20
```

Each realization writes immediately to:

```text
out/.../statistical_procedure/del120s_ac300s_1Gy_D1p0Gy/washed/
    realizationNNNN.toml
```

If needed later, obtain a nominal ensemble by replacing `--mode washed` with
`--mode nominal`.

## 3. Combine completed realizations

After the desired washed realizations have completed, calculate raw and
finite-pool-corrected range precision:

```bash
julia --project=. drivers/statistical_procedure_jobs.jl \
  --stage combine --mode washed
```

This writes:

```text
out/.../statistical_procedure/del120s_ac300s_1Gy_D1p0Gy/combined/
    washed_N100.toml
```

The combine stage checks compatible metadata and duplicate realization indices
before calculating `sigma_R`.
