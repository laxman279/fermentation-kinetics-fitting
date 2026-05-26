import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# True parameters (these are what you'll try to recover)
mu_max_true = 0.45   # 1/h
Ks_true = 0.12       # g/L
Yxs_true = 0.5       # g biomass / g substrate

# Initial conditions
X0 = 0.05   # g/L biomass
S0 = 10.0   # g/L substrate

def monod_ode(t, y, mu_max, Ks, Yxs):
    X, S = y
    mu = mu_max * S / (Ks + S) if S > 0 else 0
    dXdt = mu * X
    dSdt = -(1/Yxs) * mu * X
    return [dXdt, dSdt]

# Solve forward in time
t_span = (0, 20)
t_eval = np.linspace(0, 20, 25)
sol = solve_ivp(monod_ode, t_span, [X0, S0],
                args=(mu_max_true, Ks_true, Yxs_true),
                t_eval=t_eval, dense_output=True)

# Add Gaussian noise to simulate measurement error
np.random.seed(42)
# Use 5% relative noise per point (proportional to local value), not 5% of max
X_noisy = sol.y[0] * (1 + np.random.normal(0, 0.05, len(t_eval)))
S_noisy = sol.y[1] * (1 + np.random.normal(0, 0.05, len(t_eval)))

# Clip to non-negative
X_noisy = np.clip(X_noisy, 0, None)
S_noisy = np.clip(S_noisy, 0, None)

df = pd.DataFrame({
    'time_h': t_eval,
    'biomass_gL': X_noisy,
    'substrate_gL': S_noisy
})

df.to_csv('data/synthetic_batch_fermentation.csv', index=False)
print(df.head())
print(f"\nGround truth: mu_max={mu_max_true}, Ks={Ks_true}, Yxs={Yxs_true}")