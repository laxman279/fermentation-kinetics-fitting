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
```
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
```
## How to run

```bash
python -m venv .venv
.\.venv\Scripts\activate     # Windows PowerShell
pip install -r requirements.txt
python src/generate_data.py
jupyter notebook
```

Then open notebooks in order 01 → 04.

## Results

### Parameter recovery (Monod)

The Monod model recovered μ_max and Yxs within ~10% of ground-truth values
from synthetic batch data:

- μ_max = 0.48 h⁻¹ (ground truth: 0.45)
- Yxs   = 0.44 g/g (ground truth: 0.50)
- Ks    = 0.72 g/L (ground truth: 0.12)

Ks recovery was substantially poorer than the other parameters. This reflects
a well-known parameter identifiability limitation: when S₀ >> Ks throughout
most of the batch (here, S₀ = 10 g/L vs Ks = 0.12 g/L), the Monod term
S/(Ks+S) is approximately 1 for most of the experiment, leaving Ks weakly
constrained by the data. Reliable Ks estimation typically requires chemostat
data at varying dilution rates rather than batch experiments. The result was
robust to multiple initial guesses (the optimizer converged to the same
parameter values from different starting points), indicating this is not a
local-minimum issue but a structural property of batch data.

### Model comparison (Monod vs Haldane)

The Haldane substrate-inhibition model was fit to the same data using a
4-parameter version with Ki as the additional inhibition constant. Two
independent lines of evidence converge on the same conclusion: substrate
inhibition is not supported by the data.

1. **Ki pinned at upper bound.** The Haldane optimizer drove Ki to the upper
   bound of the search range (1000 g/L). As Ki → ∞, the Haldane term
   S²/Ki → 0 and Haldane mathematically reduces to Monod. The optimizer
   was implicitly trying to recover Monod from Haldane.

2. **Akaike Information Criterion.** ΔAIC = AIC(Haldane) − AIC(Monod) = +2.13,
   marginally favoring Monod. The Haldane fit's higher AIC primarily reflects
   the parsimony penalty (+2 per extra parameter) rather than substantially
   worse residuals. The simpler model is preferred not because it fits better,
   but because the more complex model fits no better despite an extra degree
   of freedom — the textbook signature of an unjustified parameter.

This is the expected result, since the synthetic data was generated under
pure Monod kinetics (no substrate inhibition).


## Learnings
1.Both Monod and Haldane fits achieved comparable fits to the data (final cost: Monod 0.0200, Haldane 0.0197), but with very different parameter values. This illustrates a well-known parameter identifiability problem with batch fermentation data: when S0 >> Ks throughout most of the experiment, Ks is poorly constrained, and the optimizer finds non-unique parameter combinations. Reliable Ks estimation typically requires chemostat data at varying dilution rates.


## Author

Laxman Giri — B.Tech, Biochemical Engineering and Biotechnology, IIT Delhi
