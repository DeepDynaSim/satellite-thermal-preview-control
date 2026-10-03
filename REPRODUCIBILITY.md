# Reproducibility package version 2.0.0

This package accompanies **Thermal Preview Learning and Lyapunov Safeguarding for Satellite Appendage Vibration Control**. It contains simulation evidence and editable manuscript outputs for a major revision. All datasets are computed benchmarks. Publisher PDFs and third-party book illustrations are not redistributed.

## Environment

Use Python 3.12 and the exact versions in `requirements.txt`. The numerical code uses float64, NumPy/SciPy, Numba with `fastmath=False`, and one BLAS thread. Scientific plots are exported as 420 dpi PNG and vector PDF/SVG. Independent nominal algebra/integration validation was performed in MATLAB R2024b; its report and script are included. Authoring uses python-docx, latex2mathml and Microsoft's installed MathML-to-OMML stylesheet. The supplied DOCX/PDF can be read without rerunning authoring.

From the extracted package directory, install requirements into a dedicated environment and run:

```powershell
python code/verify_package.py
python -m pip install -r requirements.txt
python code/revision_studies.py check
python code/revision_studies.py lhs --N 320 --workers 4
python code/additional_studies.py
python code/rbf_comparator.py
python code/finite_element_validation.py
python code/refinement_and_hardware.py
python code/observer_refinement.py
python code/observer_temperature_calibration.py
python code/residual_sensor_study.py
python code/timing_and_electrical.py
python code/revision_figures.py
```

Commands are ordered here for a clean sequential reproduction. The full uncertainty batch can take over an hour depending on the processor. `refinement_and_hardware.py` includes the high-mode time refinement; `high_mode_timestep.py` is an optional isolated rerun of that check. The full study creates the required report files before timing and plotting. Do not run the legacy `simulate.py --mc` command afterward unless intentionally regenerating the original 40-pair submission: it writes the original report/figure conventions.

The initial NN labels, fixed split and selected weights are in `data/learning.npz` and `data/network.npz`. They are retained exactly from the initial fully reproducible 6000-point computation; the new architecture studies use the same labels and held-out points. `simulate.py` contains the original data-generation and training code if that stage is to be regenerated in a separate folder. `data/matlab_validation.json` and `independent_validation.mat` preserve the independent original nominal check.

## Model and metric definitions

- Native baseline interval is 1/6400 s; 65-mode time checks additionally use 1/12800 and 1/25600 s.
- Revised principal comparisons integrate 120 s, with RMS including all native samples. The initial submission used a different RMS/time window and is not mixed into the revised result tables.
- LHS seed is 20261004, paired bootstrap seed 55041. Each of 320 designs is reused by LQR, analytic preview, raw NN and filtered NN.
- Every channel statistics array has rows: peak, RMS, applied saturation fraction, pre-filter clipping fraction, longest consecutive saturation duration.
- `max_nominal_true_margin` evaluates the **nominal formula at the true retained state**. It is not the exact uncertain/high-order plant derivative. Observer filtering uses estimated states; this distinction is intentional.
- Failure termination is distinct from linear instability. A saturated finite trajectory can remain below the termination threshold while the sampled linear spectrum is unstable.
- FE validation uses independently assembled Hermite elements and mass-normalized patched clamped eigenmodes with modal damping. It retains the same planar beam assumptions and is not hardware validation.
- Spatial field-aware tests require field reference information; their residual map is analytical, and the original NN learns only the mean-temperature preview component.
- CPU timings cover cached analytic preview, NN, cubic polynomial and 256-center RBF feedforward evaluation on the identified desktop CPU after warming, excluding preparation and the common feedback/observer/filter.
- Electrical accounting assumes permittivity, output resistance and loss tangent. Its spectral dielectric quantity is a finite-record proxy, not a measured power budget.

## Files and audit trail

`data/` contains per-case outputs, complete reports, design samples, training parameters, exact formula sources and verified reference metadata. `figures/` contains raster/vector figures. `code/` contains all numerical and document-building code. `Response_to_Reviewer.md` maps every critique to evidence and maps the archive deposit to its version DOI. SHA-256 manifests identify the delivered files.

## Permanent archive and licenses

Version 2.0.0 is deposited in Zenodo: https://doi.org/10.5281/zenodo.23126822.

The sole archive creator is **Caglar Uyulan**, Faculty of Engineering and Architecture, Department of Mechanical Engineering, İzmir Kâtip Çelebi University, İzmir, Türkiye. Contact: caglar.uyulan@ikcu.edu.tr.
The deposit creator metadata is distinct from the scholarly authorship and source provenance retained in the manuscript.

Generated data, figures and documentation: **CC BY 4.0**, as stated in `LICENSE_DATA.md`.
Newly authored source code: **MIT**, full terms in `LICENSE_CODE.txt`.
Source publisher PDFs, book illustrations and supplied original manuscripts are excluded from redistribution.
The actual version DOI was reserved in the authenticated Zenodo draft before file assembly. Public publication and
DOI resolution are checked at release; the local `Zenodo_Publication_Receipt.json` records that subsequent verification.
