# Fermentation Kinetics Fitting

Numerical fitting of Monod and Haldane growth kinetic models to batch microbial fermentation data, with model comparison and residual analysis.

## Motivation

Predicting and optimizing microbial cell factories requires kinetic parameters (μ_max, Ks, Ki, Yx/s) that describe how growth rate depends on substrate concentration. This project implements end-to-end pipelines for fitting these models to time-course data and evaluating goodness of fit.

## Methods

- **Models implemented:** Monod (no inhibition), Haldane (substrate inhibition), and logistic (substrate-independent) growth.
- **Numerical solver:** `scipy.integrate.solve_ivp` (Runge-Kutta 4(5), adaptive step).
- **Parameter estimation:** non-linear least squares via `scipy.optimize.minimize` with bounded parameters, with comparison against `lmfit`.
- **Model selection:** Akaike Information Criterion (AIC) and residual analysis.

## Dataset

Currently uses synthetic batch fermentation data generated from a known Monod model with added Gaussian measurement noise (5% relative). This allows validation of the fitting pipeline against ground truth. See `src/generate_data.py` and `data/README.md`.

## Repository structure
fermentation-kinetics-fitting/
├── README.md
├── requirements.txt
├── .gitignore
├── data/
│   ├── README.md
│   └── synthetic_batch_fermentation.csv
├── notebooks/
│   ├── 01_explore.ipynb
│   ├── 02_monod_fit.ipynb
│   ├── 03_haldane_fit.ipynb
│   └── 04_residual_analysis.ipynb
├── figures/
└── src/
├── generate_data.py
└── models.py

## How to run

```bash
python -m venv .venv
.\.venv\Scripts\activate     # Windows PowerShell
pip install -r requirements.txt
python src/generate_data.py
jupyter notebook
```

Then open notebooks in order 01 → 04.

## Status

- [x] Data generation
- [x] Exploratory plots
- [ ] Monod model fit
- [ ] Haldane model fit
- [ ] Residual analysis
- [ ] Final interpretation

## Author

Laxman Giri — B.Tech, Biochemical Engineering and Biotechnology, IIT Delhi