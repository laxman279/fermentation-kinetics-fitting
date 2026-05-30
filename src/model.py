"""
Kinetic models for batch fermentation.

Each model is defined as a system of ODEs returning [dX/dt, dS/dt],
where X is biomass concentration (g/L) and S is substrate concentration (g/L).
"""

import numpy as np
from scipy.integrate import solve_ivp


def monod_ode(t, y, mu_max, Ks, Yxs):
    """
    Monod growth kinetics: mu = mu_max * S / (Ks + S)

    Parameters
    ----------
    t : float
        Time (h). Required by solve_ivp signature; not used directly.
    y : array of shape (2,)
        Current state [X, S].
    mu_max : float
        Maximum specific growth rate (1/h).
    Ks : float
        Half-saturation constant (g/L).
    Yxs : float
        Biomass yield on substrate (g biomass / g substrate).

    Returns
    -------
    [dX/dt, dS/dt]
    """
    X, S = y
    S = max(S, 0)  # guard against negative substrate from solver overshoot
    mu = mu_max * S / (Ks + S) if S > 0 else 0
    dXdt = mu * X
    dSdt = -(1 / Yxs) * mu * X
    return [dXdt, dSdt] 

def simulate_monod(params, t_eval, X0, S0):
    """
    Forward-simulate the Monod model and return predicted [X, S] at t_eval.

    Returns arrays of NaN if the ODE solver fails — this allows downstream
    optimization routines to detect and avoid pathological parameter sets
    without crashing.

    Parameters
    ----------
    params : array-like of length 3
        [mu_max, Ks, Yxs]
    t_eval : array
        Time points at which to evaluate the solution (h).
    X0, S0 : float
        Initial biomass and substrate concentrations (g/L).

    Returns
    -------
    X_pred, S_pred : arrays of length len(t_eval)
    """
    mu_max, Ks, Yxs = params
    sol = solve_ivp(
        monod_ode,
        t_span=(t_eval.min(), t_eval.max()),
        y0=[X0, S0],
        args=(mu_max, Ks, Yxs),
        t_eval=t_eval,
        method='RK45',
        rtol=1e-6,
        atol=1e-9,
    )
    if not sol.success:
        # Return arrays of NaN so the fitting routine knows this is a bad parameter set
        return np.full_like(t_eval, np.nan), np.full_like(t_eval, np.nan)
    return sol.y[0], sol.y[1]


def residuals_monod(params, t_data, X_data, S_data, X0, S0):
    """
    Compute residuals between predicted and observed [X, S] for given parameters.

    Returns a flat array of length 2*N where N = len(t_data),
    suitable for scipy.optimize.least_squares.
    """
    X_pred, S_pred = simulate_monod(params, t_data, X0, S0)

    # Normalize each variable by its scale so X and S contribute comparably
    X_scale = np.nanmax(X_data) if np.nanmax(X_data) > 0 else 1
    S_scale = np.nanmax(S_data) if np.nanmax(S_data) > 0 else 1

    X_res = (X_pred - X_data) / X_scale
    S_res = (S_pred - S_data) / S_scale

    return np.concatenate([X_res, S_res])

def haldane_ode(t, y, mu_max, Ks, Ki, Yxs):
    """
    Haldane growth kinetics with substrate inhibition:
        mu = mu_max * S / (Ks + S + S^2 / Ki)

    The S^2/Ki term penalizes growth at high substrate concentrations.
    As Ki -> infinity, Haldane reduces to Monod.

    Parameters
    ----------
    t : float (required by solve_ivp signature)
    y : [X, S]
    mu_max : float, max specific growth rate (1/h)
    Ks : float, half-saturation constant (g/L)
    Ki : float, substrate inhibition constant (g/L). Larger => less inhibition.
    Yxs : float, biomass yield (g/g)
    """
    X, S = y
    S = max(S, 0)
    if S > 0:
        mu = mu_max * S / (Ks + S + S**2 / Ki)
    else:
        mu = 0
    dXdt = mu * X
    dSdt = -(1 / Yxs) * mu * X
    return [dXdt, dSdt]


def simulate_haldane(params, t_eval, X0, S0):
    """Forward-simulate the Haldane model. params = [mu_max, Ks, Ki, Yxs]"""
    import numpy as np
    from scipy.integrate import solve_ivp

    mu_max, Ks, Ki, Yxs = params
    sol = solve_ivp(
        haldane_ode,
        t_span=(t_eval.min(), t_eval.max()),
        y0=[X0, S0],
        args=(mu_max, Ks, Ki, Yxs),
        t_eval=t_eval,
        method='RK45',
        rtol=1e-6,
        atol=1e-9,
    )
    if not sol.success:
        return np.full_like(t_eval, np.nan), np.full_like(t_eval, np.nan)
    return sol.y[0], sol.y[1]


def residuals_haldane(params, t_data, X_data, S_data, X0, S0):
    """Residual function for Haldane, scaled by max of each variable."""
    import numpy as np
    X_pred, S_pred = simulate_haldane(params, t_data, X0, S0)
    X_scale = np.nanmax(X_data) if np.nanmax(X_data) > 0 else 1
    S_scale = np.nanmax(S_data) if np.nanmax(S_data) > 0 else 1
    X_res = (X_pred - X_data) / X_scale
    S_res = (S_pred - S_data) / S_scale
    return np.concatenate([X_res, S_res])