# Vendored reference docs

Snapshots of the upstream contracts this repo builds against, copied here so CryspBrainSim reads as
a self-contained whole. The upstream files remain the source of truth; refresh these copies when the
upstream repos advance.

Captured 2026-07-05; `PRODUCTS.md` and `SCHEMA.md` refreshed 2026-07-10 for the two-scanner
production (PTCryspMC.jl @ `27aab25`; published copies at the PtCryspProds root are the source).

| Local file | Source repo @ SHA | Upstream path |
|---|---|---|
| `PRODUCTS.md` | PTCryspMC.jl @ `27aab25` (main) | published `PtCryspProds/README.md` |
| `data_generation_strategy.md` | PTCryspMC.jl @ `23b1674` (main) | `dev/data_generation_strategy.md` |
| `range_verification_recipe.md` | PTCryspMC.jl @ `23b1674` (main) | `docs/range_verification_recipe.md` |
| `SCHEMA.md` | PTCryspMC.jl @ `27aab25` (main) | published `PtCryspProds/SCHEMA.md` |
| `RecoCryspUse.md` | RecoCrysp @ `a6d900e` (main) | `statusmd/RecoCryspUse.md` |

To refresh: re-copy from the upstream paths above and update the SHAs and date in this table.
