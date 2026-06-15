# Dataset

`synthetic_batch_fermentation.csv` — synthetic time-course data generated from a Monod growth model using `src/generate_data.py`.

## Ground-truth parameters used to generate the data

- μ_max = 0.45 h⁻¹
- Ks = 0.12 g/L
- Yx/s = 0.5 g biomass / g substrate
- Initial biomass X₀ = 0.05 g/L
- Initial substrate S₀ = 10.0 g/L
- Measurement noise: 5% relative Gaussian(proportional to local value)

## Columns

- `time_h` — time in hours
- `biomass_gL` — biomass concentration in g/L
- `substrate_gL` — substrate concentration in g/L
