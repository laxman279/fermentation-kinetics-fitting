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

Uses synthetic batch fermentation data generated from a known Monod model with added Gaussian measurement noise (5% relative). This allows validation of the fitting pipeline against ground truth. See `src/generate_data.py` and `data/README.md`.

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

## Parameter rationale

The ground-truth parameters approximate a *Saccharomyces cerevisiae* aerobic
batch fermentation on glucose:

- μ_max = 0.45 1/h corresponds to a doubling time of ~1.5 h, typical for
  yeast on glucose under aerobic conditions.
- Ks = 0.12 g/L is in the standard range for yeast on glucose (literature: 0.1–1 g/L).
- Yxs = 0.50 g/g reflects aerobic respiratory growth (theoretical maximum ~0.55 g/g).
- S0 = 10 g/L (1% glucose) is a standard laboratory fermentation substrate
  concentration.
- X0 = 0.05 g/L corresponds to a typical inoculum OD600 ≈ 0.1.

The ratio S0/Ks ≈ 83 is realistic for industrial fermentations but is also
what creates the Ks identifiability challenge demonstrated in this project.

## Results

### Synthetic-data ground truth

The dataset was generated from a known Monod model with relative Gaussian noise:

```text
mu_max_true = 0.45 1/h
Ks_true     = 0.12 g/L
Yxs_true    = 0.50 g/g
X0          = 0.05 g/L
S0          = 10.0 g/L
noise       = 5% relative Gaussian
n_data      ≈ 31 time points
```

Because the ground truth is known, the fitted parameters can be evaluated directly rather than only judged by visual fit quality.

### Parameter recovery: Monod model

The Monod model recovered μ_max and Yxs reasonably well, but not Ks:

```text
mu_max = 0.4802 1/h   (truth 0.45 → within ~7%)
Ks     = 0.7245 g/L   (truth 0.12 → off by ~6×)
Yxs    = 0.4441 g/g   (truth 0.50 → within ~12%)
final cost (SSR/2) = 0.060898
```

The key result is not that the fit failed. The overall fit is clean, but Ks is poorly identified. Four observations support this interpretation:

1. Ks is far from the ground truth despite a good overall fit.
2. The result is robust to initial guess; different starting points converge to the same value.
3. The problem is model-independent; Monod and Haldane return nearly identical Ks values.
4. Substrate residuals show a systematic negative bias, consistent with overestimated Ks and slower predicted substrate consumption.

The root cause is that S0 = 10 g/L is much larger than Ks = 0.12 g/L. For most of the batch, the Monod term S/(Ks + S) is close to 1, so Ks has little influence on the trajectory. Reliable Ks estimation would require experiments that probe low-substrate regimes, such as chemostat experiments at varying dilution rates.

### Model comparison: Monod vs Haldane

The Haldane model gave almost the same fit, but its extra inhibition parameter was not supported:

```text
mu_max = 0.4863 1/h
Ks     = 0.7610 g/L
Ki     = 1000.0 g/L   (pinned at upper bound)
Yxs    = 0.4441 g/g
final cost (SSR/2) = 0.061030
```

As Ki approaches infinity, the Haldane term S²/Ki approaches zero, so the Haldane model mathematically collapses back to Monod. The optimizer driving Ki to the upper bound means it was trying to remove substrate inhibition from the model.

The AIC comparison also favors Monod:

```text
AIC (Monod, 3 params)   = -380.42
AIC (Haldane, 4 params) = -378.28
ΔAIC (Haldane - Monod)  = +2.13
```

The higher AIC for Haldane mainly reflects the penalty for the extra parameter. Haldane also has a slightly worse final cost than Monod (0.061030 vs 0.060898), so the extra parameter does not improve the fit. This is expected because the synthetic data was generated under pure Monod kinetics.

### Residual diagnostics

The residual statistics were:

```text
Monod biomass:      mean = -0.0183, std = 0.1684   (|mean|/std ≈ 0.11)
Monod substrate:    mean = -0.3327, std = 0.4917   (|mean|/std ≈ 0.68)
Haldane biomass:    mean = -0.0181, std = 0.1685
Haldane substrate:  mean = -0.3335, std = 0.4921
```

Biomass residuals are centered near zero and show no strong time-dependent trend. Substrate residuals have a stronger negative bias, consistent with the Ks identifiability issue.

Shapiro-Wilk normality tests gave:

```text
Monod biomass:      W = 0.9368, p = 0.0673  → fail to reject normality
Monod substrate:    W = 0.7526, p < 0.0001  → reject normality
Haldane biomass:    W = 0.9373, p = 0.0696  → fail to reject normality
Haldane substrate:  W = 0.7540, p < 0.0001  → reject normality
```

The biomass residuals are consistent with the Gaussian-noise assumption used in ordinary least squares. The substrate residuals are non-normal because the data were generated with 5% relative noise. This makes the residual variance depend on substrate concentration and shrink near depletion, violating the constant-variance assumption. The point estimates remain useful, but uncertainty estimates and AIC would be more rigorous under weighted least squares, with residuals weighted by inverse expected variance.

### Bootstrap identifiability analysis

Bootstrap resampling showed that Ks is much less stable than Yxs:

```text
Parameter    95% CI width
mu_max       0.6259
Ks           6.0435
Yxs          0.0184
```

Ks and μ_max were strongly coupled across bootstrap resamples (r = +0.96). This explains the identifiability problem: the optimizer can move along a broad μ_max-Ks ridge while producing similar batch trajectories. Yxs was nearly uncorrelated with both parameters (r ≈ 0.1), indicating that it is independently identifiable from this dataset.

## Learnings

1. Batch fermentation data can recover μ_max and Yxs reasonably well, but Ks is weakly constrained when S0 >> Ks.
2. The poor Ks estimate is not a local-minimum problem. It is structural: the data do not contain enough low-substrate information to identify Ks.
3. Haldane is not justified for this dataset. Ki is pinned at the upper bound, the final cost is slightly worse than Monod, and AIC favors the simpler Monod model.
4. Residual diagnostics reveal that biomass fits satisfy the Gaussian-noise assumption better than substrate fits.
5. Weighted least squares would be a more statistically rigorous next step because the dataset uses relative measurement noise.
6. Bootstrap analysis confirms the mechanism of non-identifiability: μ_max and Ks are highly correlated, while Yxs remains independently identifiable.

## Written by-

Laxman Giri — B.Tech, Biochemical Engineering and Biotechnology, IIT Delhi
