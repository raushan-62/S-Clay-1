# S-CLAY1 Material Point Workbench

A local, browser-based triaxial material-point simulator for the simplified
S-CLAY1 model by Wheeler, Näätänen, Karstunen, and Lojander (2003). The web
interface can also overlay a Modified Cam Clay (MCC) baseline for the same
strain path.

## Files

- `index.html` — the workbench interface.
- `styles.css` — responsive page styling.
- `app.js` — browser-side integration and charts. It has no runtime library
  dependencies and sends no inputs or results to a server.
- `s_clay1.py` — Python material-point implementation of the same simplified
  S-CLAY1 formulation.

## Start the website

1. Keep `index.html`, `styles.css`, and `app.js` together in one folder.
2. Open `index.html` in a current desktop or mobile browser.
3. Enter soil constants and an admissible initial state, or use the defaults.
4. Choose a strain path and click **Run simulation**.
5. Use **Export CSV** to save the point-by-point history or the chart button to
   save a PNG of the plotted response.

The website is static: it does not need Python, a package install, or a
backend. To share the app as a hosted link, publish all three website files
(`index.html`, `styles.css`, and `app.js`) together as a static site. Sharing
only the HTML file will omit the styling and solver script.

For GitHub Pages, place those files at the root of a GitHub repository, then
open **Settings → Pages**, choose deployment from the `main` branch and `/`
(root), and save. GitHub will provide the website URL after the deployment
finishes. A public Pages site makes the source files visible to anyone with
access to the repository.

## Model and sign conventions

- `p'` and `q` are effective stresses in kPa; compression is positive.
- `q` is signed: negative `q` represents triaxial extension.
- `eps_v` and `eps_q` are compression-positive, work-conjugate strains. The
  work increment is `p' d(eps_v) + q d(eps_q)`.
- `v` is specific volume, `1 + e`.
- Yield function:

  ```text
  f = (q - alpha p')^2 - (M^2 - alpha^2)(p'_m - p')p' = 0
  ```

- The model uses associated flow, isotropic elastic response, yield-size
  hardening, and rotational hardening with `chi_v = 3 eta / 4` and
  `chi_d = eta / 3`, where `eta = q / p'`.
- Integration uses an elastic predictor, bisection to locate yield, an explicit
  plastic multiplier, and a size-only drift correction. The maximum normalized
  yield drift before correction is shown and recorded.

This is the simplified triaxial-invariant version, not the tensor model for
general 3D stress paths. A single value of `M` is used in compression and
extension.

## Inputs

### Soil constants

| Input | Meaning |
|---|---|
| `lambda` (`λ`) | Virgin compression slope |
| `kappa` (`κ`) | Elastic swelling slope; must be positive and less than `λ` |
| `nu` (`ν'`) | Elastic Poisson ratio |
| `M` | Critical-state stress ratio |
| `beta` (`β`) | Relative shear-strain contribution to rotational hardening |
| `mu` (`μ`) | Rotational hardening rate |

### Initial state

| Input | Meaning |
|---|---|
| `p'0`, `q0` | Initial effective stress state, kPa |
| `p'm0` | Initial yield-surface size, kPa |
| `v0` | Initial specific volume; must exceed 1 |
| `alpha0` (`α0`) | Initial yield-surface inclination |
| `etaK0` | Stress ratio used by the **Set α** estimate |

**Set α** calculates the normally consolidated K0 estimate
`alpha_K0 = eta_K0 - 3/4 (M^2 - eta_K0^2)`. **OCR 1** sets `p'm0` so that the
current stress lies on the configured initial yield surface. Check the
admissibility indicator below the fields before running.

The displayed defaults are demonstration values, not a calibrated soil data
set. Replace them with the values reported for the soil and paper used in your
project.

## Strain paths

All paths apply prescribed increments to the material point.

| Protocol | Strain path |
|---|---|
| Undrained triaxial compression | `eps_v = 0`; positive total `eps_q` |
| Undrained triaxial extension | `eps_v = 0`; negative total `eps_q` |
| Cyclic undrained triaxial | Triangular cycles `0 -> +amplitude -> -amplitude -> 0`; choose cycles and steps per leg |
| Isotropic compression/swelling | `eps_q = 0`; enter total `eps_v` (positive is compression) |
| K0 / oedometer | `d eps_q = 2 d eps_v / 3`, corresponding to `d eps_r = 0` |
| Custom linear path | Enter total `eps_v` and `eps_q` for a fixed-ratio path |

**Load steps** control the prescribed path resolution. **Internal step limit**
sets the maximum strain magnitude in an integration substep. Smaller internal
steps usually reduce explicit-integration drift, at the cost of run time.

## Modified Cam Clay comparison

The **Overlay Modified Cam Clay** option is on by default. The baseline keeps
`lambda`, `kappa`, `nu`, `M`, initial stresses, and specific volume the same,
and sets `alpha = 0`, `beta = 0`, and `mu = 0`. This removes the anisotropic
yield inclination and rotational hardening, giving the MCC special case in
this implementation. Both model histories are plotted together and included
in the CSV.

If the chosen MCC yield size would place the starting stress outside its
surface, the app raises only the MCC baseline's initial `p'm` to the minimum
size that contains that starting stress. The comparison note reports this
adjustment. S-CLAY1 inputs are never changed by this adjustment.

## Results and exports

The charts show:

1. Deviator stress `q` against deviatoric strain `eps_q`.
2. Effective stress path `q` against `p'`, with initial and final yield
   surfaces.
3. Specific volume `v` against `p'`.
4. Yield inclination `alpha` against `eps_q`.
5. For undrained paths, estimated excess pore pressure against `eps_q`.

The excess pore pressure estimate assumes constant total cell pressure:
`delta_u = p'0 - p' + (q - q0) / 3`. It is an estimate under that boundary
condition; it is not suitable for every undrained stress-control setup.

The CSV contains step number and path fraction, stresses, strains, specific
volume, yield size, inclination, yield function, pre-correction drift,
cumulative plastic multiplier, substep count, and estimated excess pore
pressure. When MCC comparison is enabled, MCC stress, volume, yield size, and
yield-function columns are added alongside the S-CLAY1 history.

## Python interface

The Python routine uses only the standard library. Example:

```python
from s_clay1 import Parameters, initial_state, stress_update_substepped

soil = Parameters(lam=0.20, kappa=0.04, nu=0.30,
                  M=1.20, beta=0.70, mu=50.0)
state = initial_state(v=2.0, p_m=220.0, alpha=0.42)
sigma = (100.0, 0.0)  # (p', q), kPa
strain = (0.0, 0.0)   # (eps_v, eps_q)

sigma, strain, state = stress_update_substepped(
    sigma, strain, (0.0, 0.001), state, soil
)
```

Available functions:

- `stress_update(sigma, eps, d_eps, state, params)` — one strain increment.
- `stress_update_substepped(...)` — one increment split into smaller steps.
- `run_strain_path(...)` — a sequence of prescribed increments.
- `yield_function(...)` — evaluate the S-CLAY1 yield function.
- `alpha_k0(M, eta_k0)` — estimate initial inclination for K0 loading.

The return value from each update is `(sigma_new, eps_new, state_new)`.

## Important limits

- The source paper cautions that this simple isotropic-elastic S-CLAY1 model is
  unsuitable for cyclic loading. The cyclic option is included for exploratory
  course work and should not be presented as a validated cyclic prediction.
- The MCC overlay is a comparison baseline, not an independently calibrated
  MCC implementation or a separate constitutive integration scheme.
- A successful run only shows that the selected inputs and path integrated;
  it does not validate the model against measured soil behavior.
- The course verification sweep, source-paper parameter reproduction,
  digitized-data validation, and engineering interpretation remain separate
  project tasks.

## Troubleshooting

- **Initial stress outside yield:** increase `p'm0`, or use **OCR 1** to place
  the surface through the current stress; confirm `|alpha| < M`.
- **`lambda` / `kappa` error:** use `lambda > kappa > 0`.
- **Negative pressure or unstable update:** reduce the prescribed load-step
  size by increasing the number of path steps, and/or reduce the internal step
  limit.
- **MCC line absent:** ensure the overlay checkbox is selected. If the MCC
  update fails for the selected path, the app keeps the S-CLAY1 result and
  reports the baseline error.
- **Website looks unstyled or Run does nothing:** make sure all three web files
  are together in one folder and reopen `index.html`.

## Reference

Wheeler, S. J., Näätänen, A., Karstunen, M., and Lojander, M. (2003).
“An anisotropic elastoplastic model for soft clays.” *Canadian Geotechnical
Journal*, 40, 403–418. DOI: 10.1139/T02-119.
