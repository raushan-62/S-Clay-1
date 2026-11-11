"""Triaxial-invariant S-CLAY1 constitutive model (Wheeler et al., 2003).

Stress convention: compression-positive effective stresses ``(p, q)`` in kPa.
Strains: compression-positive work-conjugate invariants ``(eps_v, eps_q)``;
plastic work is ``p*d_eps_v + q*d_eps_q``.  ``q`` is signed, so extension is
negative.  This is the simplified triaxial model, not the tensor generalisation.

The update is explicit and evaluates moduli/gradients at each substep start.
It uses an elastic predictor, bisection to locate first yield, a plastic
multiplier for the remainder, and a size-only drift correction.  The maximum
pre-correction normalized drift is accumulated in State.max_normalized_drift.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from math import ceil, exp, isfinite
from typing import Iterable, Tuple


Vector = Tuple[float, float]


@dataclass(frozen=True)
class Parameters:
    lam: float                 # virgin compression slope lambda
    kappa: float               # swelling slope kappa
    nu: float                  # elastic Poisson ratio
    M: float                   # triaxial critical-state stress ratio (positive)
    beta: float                # relative shear-strain contribution to rotation
    mu: float                  # rotational hardening rate
    p_floor: float = 1.0e-6
    alpha_limit_ratio: float = 0.999
    yield_tolerance: float = 1.0e-10
    bisection_iterations: int = 60

    def validate(self) -> None:
        if not (self.lam > self.kappa > 0):
            raise ValueError("Require lambda > kappa > 0.")
        if not (-0.99 < self.nu < 0.49):
            raise ValueError("Poisson ratio must lie between -0.99 and 0.49.")
        if self.M <= 0 or self.beta < 0 or self.mu < 0:
            raise ValueError("Require M > 0, beta >= 0, and mu >= 0.")
        if self.p_floor <= 0 or not (0 < self.alpha_limit_ratio < 1):
            raise ValueError("Invalid pressure floor or alpha limit ratio.")


@dataclass(frozen=True)
class State:
    v: float                    # specific volume
    p_m: float                  # yield-surface size p'_m, kPa
    alpha: float                # yield-surface inclination
    max_normalized_drift: float = 0.0
    plastic_multiplier_sum: float = 0.0
    substeps: int = 0
    eps_v_plastic: float = 0.0  # accumulated associated plastic volumetric strain
    eps_q_plastic: float = 0.0  # accumulated associated plastic deviatoric strain


def initial_state(v: float, p_m: float, alpha: float = 0.0) -> State:
    if v <= 1.0 or p_m <= 0:
        raise ValueError("Require v > 1 and p_m > 0.")
    return State(v=v, p_m=p_m, alpha=alpha)


def alpha_k0(M: float, eta_k0: float) -> float:
    """S-CLAY1 normally-consolidated inclination for one-dimensional loading."""
    if M <= abs(eta_k0):
        raise ValueError("K0 stress ratio must satisfy |eta_k0| < M.")
    return eta_k0 - 0.75 * (M * M - eta_k0 * eta_k0)


def yield_function(p: float, q: float, p_m: float, alpha: float, M: float) -> float:
    """Sheared-ellipse yield function f=(q-alpha*p)^2-(M^2-alpha^2)(p_m-p)p."""
    return (q - alpha * p) ** 2 - (M * M - alpha * alpha) * (p_m - p) * p


def _yield_gradient(p: float, q: float, p_m: float, alpha: float, M: float):
    x = q - alpha * p
    fp = -2.0 * alpha * x - (M * M - alpha * alpha) * (p_m - 2.0 * p)
    fq = 2.0 * x
    fpm = -(M * M - alpha * alpha) * p
    fa = -2.0 * p * x + 2.0 * alpha * (p_m - p) * p
    return fp, fq, fpm, fa


def _moduli(p: float, v: float, prm: Parameters):
    K = v * p / prm.kappa
    G = 3.0 * K * (1.0 - 2.0 * prm.nu) / (2.0 * (1.0 + prm.nu))
    return K, 3.0 * G


def _elastic(p: float, q: float, v: float, deps: Vector, prm: Parameters):
    K, threeG = _moduli(p, v, prm)
    p1 = p + K * deps[0]
    q1 = q + threeG * deps[1]
    return p1, q1


def _valid(p: float, q: float, s: State, prm: Parameters) -> None:
    if not all(isfinite(x) for x in (p, q, s.v, s.p_m, s.alpha)):
        raise FloatingPointError("Non-finite state produced.")
    if p <= prm.p_floor:
        raise FloatingPointError(f"Mean effective stress reached p={p:g} kPa.")
    if s.v <= 1.0 or s.p_m <= prm.p_floor:
        raise FloatingPointError("Specific volume or yield size left its admissible range.")
    if abs(s.alpha) >= prm.alpha_limit_ratio * prm.M:
        raise FloatingPointError("Yield inclination approached |alpha| = M.")


def _one_substep(p: float, q: float, eps: Vector, deps: Vector, state: State,
                 prm: Parameters, *, locate_yield: bool = True):
    """One explicit stress-point increment; return p,q,eps,state and drift."""
    v0, pm0, a0 = state.v, state.p_m, state.alpha
    f0 = yield_function(p, q, pm0, a0, prm.M)
    scale0 = max(p * p, 1.0)
    if f0 > prm.yield_tolerance * scale0:
        raise ValueError("Initial stress lies outside the current yield surface.")

    p_trial, q_trial = _elastic(p, q, v0, deps, prm)
    _valid(p_trial, q_trial, state, prm)
    f_trial = yield_function(p_trial, q_trial, pm0, a0, prm.M)
    # Entirely elastic step, or elastic unloading into the surface.
    if not locate_yield or f_trial <= prm.yield_tolerance * max(p_trial*p_trial, 1.0):
        s1 = replace(state, v=v0 * exp(-deps[0]))
        _valid(p_trial, q_trial, s1, prm)
        return p_trial, q_trial, (eps[0] + deps[0], eps[1] + deps[1]), s1, 0.0

    # Locate first intersection along the prescribed strain increment.
    lo, hi = 0.0, 1.0
    for _ in range(prm.bisection_iterations):
        mid = 0.5 * (lo + hi)
        pp, qq = _elastic(p, q, v0, (deps[0]*mid, deps[1]*mid), prm)
        fm = yield_function(pp, qq, pm0, a0, prm.M)
        if fm > 0.0:
            hi = mid
        else:
            lo = mid
    fraction = 0.5 * (lo + hi)
    pe, qe = _elastic(p, q, v0, (deps[0]*fraction, deps[1]*fraction), prm)
    plastic_deps = (deps[0]*(1.0-fraction), deps[1]*(1.0-fraction))

    fp, fq, fpm, fa = _yield_gradient(pe, qe, pm0, a0, prm.M)
    K, threeG = _moduli(pe, v0, prm)
    # Internal-variable rates per unit plastic multiplier.
    eta = qe / pe
    alpha_rate = prm.mu * (
        (0.75*eta - a0) * max(fp, 0.0)
        + prm.beta * (eta/3.0 - a0) * abs(fq)
    )
    pm_rate = pm0 * v0 / (prm.lam - prm.kappa) * fp
    H = -(fpm * pm_rate + fa * alpha_rate)
    denominator = fp*fp*K + fq*fq*threeG + H
    numerator = fp*K*plastic_deps[0] + fq*threeG*plastic_deps[1]
    if not isfinite(denominator) or denominator <= 0:
        raise FloatingPointError("Non-positive or non-finite plastic consistency denominator.")
    dlam = numerator / denominator
    if dlam < -prm.yield_tolerance:
        # A negative multiplier indicates unloading from the newly reached point.
        s1 = replace(state, v=v0*exp(-deps[0]))
        _valid(p_trial, q_trial, s1, prm)
        return p_trial, q_trial, (eps[0]+deps[0], eps[1]+deps[1]), s1, 0.0
    dlam = max(0.0, dlam)

    # Explicit stress, size and rotation updates using start-of-plastic-step data.
    p1 = pe + K*(plastic_deps[0] - dlam*fp)
    q1 = qe + threeG*(plastic_deps[1] - dlam*fq)
    evp = dlam*fp
    edp = dlam*fq
    v1 = v0 * exp(-deps[0])
    pm1 = pm0 + dlam*pm_rate
    a1 = a0 + dlam*alpha_rate
    # Keep the ellipse non-degenerate in finite precision.
    abound = prm.alpha_limit_ratio * prm.M
    a1 = min(abound, max(-abound, a1))
    s_uncorrected = replace(state, v=v1, p_m=pm1, alpha=a1,
                            plastic_multiplier_sum=state.plastic_multiplier_sum+dlam,
                            eps_v_plastic=state.eps_v_plastic+evp,
                            eps_q_plastic=state.eps_q_plastic+edp)
    _valid(p1, q1, s_uncorrected, prm)

    # Record normalized drift before correction, using the initial pressure scale.
    f_before = yield_function(p1, q1, pm1, a1, prm.M)
    drift = abs(f_before) / max(p*p, 1.0)
    # Correct size only. Holding the computed alpha fixed preserves the rotation
    # calculated from the explicit plastic increment while returning to f=0.
    denom = (prm.M*prm.M-a1*a1)*p1
    if denom <= 0:
        raise FloatingPointError("Cannot apply size drift correction for degenerate yield ellipse.")
    pm_corrected = p1 + (q1-a1*p1)**2 / denom
    s1 = replace(s_uncorrected, p_m=pm_corrected,
                 max_normalized_drift=max(state.max_normalized_drift, drift))
    _valid(p1, q1, s1, prm)
    return p1, q1, (eps[0]+deps[0], eps[1]+deps[1]), s1, drift


def stress_update(sigma: Vector, eps: Vector, d_eps: Vector, state: State,
                  params: Parameters):
    """Update one strain increment.

    Args:
        sigma: (p_effective, q), compression positive, kPa.
        eps: accumulated (eps_v, eps_q), compression positive.
        d_eps: imposed (d_eps_v, d_eps_q).
        state: ``State(v, p_m, alpha, ...)``.
        params: S-CLAY1 soil constants.
    Returns:
        ``(sigma_new, eps_new, state_new)`` with the project brief interface.
    """
    params.validate()
    p, q = map(float, sigma)
    eps = tuple(map(float, eps))
    d_eps = tuple(map(float, d_eps))
    if len(eps) != 2 or len(d_eps) != 2:
        raise ValueError("sigma, eps, and d_eps must each have two components.")
    _valid(p, q, state, params)
    p1, q1, eps1, state1, _ = _one_substep(p, q, eps, d_eps, state, params)
    return (p1, q1), eps1, replace(state1, substeps=state1.substeps+1)


def stress_update_substepped(sigma: Vector, eps: Vector, d_eps: Vector, state: State,
                             params: Parameters, max_strain_per_substep: float = 2.5e-4,
                             max_substeps: int = 10000):
    """Apply a strain increment in equal substeps for improved explicit stability."""
    if max_strain_per_substep <= 0 or max_substeps < 1:
        raise ValueError("Invalid substep controls.")
    n = max(1, ceil(max(abs(d_eps[0]), abs(d_eps[1])) / max_strain_per_substep))
    if n > max_substeps:
        raise RuntimeError(f"Required {n} substeps; maximum is {max_substeps}.")
    de = (d_eps[0]/n, d_eps[1]/n)
    p, q = sigma
    e = eps
    s = state
    for _ in range(n):
        p, q, e, s, _ = _one_substep(p, q, e, de, s, params)
    return (p, q), e, replace(s, substeps=s.substeps+n)


def run_strain_path(sigma0: Vector, eps0: Vector, state0: State,
                    increments: Iterable[Vector], params: Parameters,
                    max_strain_per_substep: float = 2.5e-4,
                    max_substeps: int = 10000):
    """Run a sequence of imposed strain increments; return all material-point states."""
    records = [(sigma0, eps0, state0)]
    sigma, eps, state = sigma0, eps0, state0
    for de in increments:
        sigma, eps, state = stress_update_substepped(
            sigma, eps, de, state, params, max_strain_per_substep, max_substeps
        )
        records.append((sigma, eps, state))
    return records


__all__ = ["Parameters", "State", "initial_state", "alpha_k0", "yield_function",
           "stress_update", "stress_update_substepped", "run_strain_path"]

