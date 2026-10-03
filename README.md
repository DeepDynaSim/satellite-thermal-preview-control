# Satellite thermal preview control

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23126822.svg)](https://doi.org/10.5281/zenodo.23126822)
[![Package verification](https://github.com/DeepDynaSim/satellite-thermal-preview-control/actions/workflows/verify-package.yml/badge.svg)](https://github.com/DeepDynaSim/satellite-thermal-preview-control/actions/workflows/verify-package.yml)

Reproducible simulation code and the final English manuscript for **Thermal Preview Learning and Lyapunov Safeguarding for Satellite Appendage Vibration Control**.
The method combines moving thermal-equilibrium coordinates, classical finite-horizon LQ preview, nonlinear thermal forecasting,
neural feedforward compression, a conditional nominal Lyapunov filter and residual strain-rate damping.

**Repository and archive creator:** Caglar Uyulan, Faculty of Engineering and Architecture, Department of Mechanical Engineering,
İzmir Kâtip Çelebi University, İzmir, Türkiye. Contact: [caglar.uyulan@ikcu.edu.tr](mailto:caglar.uyulan@ikcu.edu.tr).
The manuscript's scholarly authorship and cited source provenance are retained in the manuscript.

## Final manuscript and permanent archive

- [Final PDF manuscript](Thermal_Preview_Learning_and_Lyapunov_Safeguarding_Revised.pdf)
- [Editable Word manuscript](Thermal_Preview_Learning_and_Lyapunov_Safeguarding_Revised.docx)
- [Permanent Zenodo version 2.0.0](https://doi.org/10.5281/zenodo.23126822)
- [Point-by-point reviewer response](Response_to_Reviewer.md) and [verified reference DOI ledger](DOI_Verification_Ledger.md)
- [Public archive verification](Zenodo_Publication_Verification.md)

The final manuscript has 30 pages, 59 numbered display equations, 144 editable Office Math objects and 11 figures.
The frozen 41.9 MB reproducibility ZIP is available in the permanent Zenodo record and the GitHub v2.0.0 release.

## Evidence and interpretation

All datasets are numerical benchmarks. Five elastic modes are used for synthesis and 37 for the principal evaluation.
The nominal filtered neural preview reduces the common 120 s cost by 56.7% at eclipse exit and 56.8% at entry relative to thermal LQR.
Across 320 paired Latin-hypercube cases, the median filtered-NN/LQR cost ratio is 0.8253, with paired-bootstrap 95% interval [0.8067, 0.8391].
The paper also includes horizon/component/architecture ablations, analytical, cubic-polynomial and 256-center RBF comparators,
independent Hermite finite elements, spatial two-face heat calculations, observer/noise/delay tests, warmed CPU timings and illustrative electrical accounting.

Exact optimality applies to the specified nominal unconstrained finite-dimensional preview problem.
The Lyapunov filter enforces a conditional nominal dissipation inequality; it does not certify arbitrary unmodeled spacecraft dynamics.
Failed sampled high-mode, observer and residual-delay configurations are retained. There is no hardware or flight validation.

## Quick integrity check

```bash
git clone https://github.com/DeepDynaSim/satellite-thermal-preview-control.git
cd satellite-thermal-preview-control
python code/verify_package.py
```

This standard-library check verifies every file in `SHA256SUMS.txt`, the 320 four-controller paired cases,
20 verified reference DOI records, 59 display equations, 144 native Word math objects and the complete visual-QA record.
It does not rerun the full simulations. GitHub Actions executes the same integrity check.

## Numerical reproduction

Use Python 3.12. The exact scientific/authoring dependency versions are pinned in `requirements.txt`.

```bash
python -m venv .venv
# Activate .venv/bin/activate on Linux/macOS, or .venv\Scripts\Activate.ps1 on Windows.
python -m pip install -r requirements.txt
python code/revision_studies.py check
```

Run the complete, ordered study commands in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).
The uncertainty batch can take over an hour depending on the processor. Numerical reruns intentionally regenerate report files;
verify the frozen downloaded repository before rerunning. The original training labels, split and selected NN weights are included.
Float64 arithmetic, `fastmath=False`, deterministic seeds and one BLAS thread are used.

MATLAB R2024b independently checked the nominal CARE, HJB identity and trajectory integration; see
[independent_validate.m](code/independent_validate.m) and [matlab_validation.json](data/matlab_validation.json).
Document regeneration additionally requires Microsoft's installed MathML-to-OMML stylesheet and Word for the native export;
reading the delivered Word/PDF does not require that authoring environment.

## Repository map

| Location | Contents |
| --- | --- |
| `code/` | Python numerical/plot/authoring scripts and independent MATLAB validation |
| `data/` | Training weights, per-case traces, complete computed reports and audit records |
| `data/lhs_cases/` | All 320 paired uncertainty-case JSON records |
| `data/references/` | Verified canonical Crossref DOI records |
| `figures/` | 11 original scientific figures in 420 dpi PNG and vector PDF/SVG |
| Root Word/PDF files | Final manuscript including the real archive DOI |
| `SHA256SUMS.txt` | GitHub repository payload integrity manifest |
| `ZENODO_SHA256SUMS.txt` | Manifest of the frozen published Zenodo package |
| `Archive_SHA256.txt` | SHA-256 of the unchanged Zenodo reproducibility ZIP |

The computed reports, figures, numerical source and manuscript match the published Zenodo payload.
The repository adds navigation, citation metadata, contributor instructions and CI; its own manifest therefore differs from the frozen ZIP manifest.
Historical audit records describe the time they were assembled; `Zenodo_Publication_Receipt.json` records the subsequent verified public release.

## Licenses and citation

Newly authored source code: **MIT**, [LICENSE](LICENSE) / [LICENSE_CODE.txt](LICENSE_CODE.txt).
Generated data, figures and documentation: **CC BY 4.0**, [LICENSE_DATA.md](LICENSE_DATA.md).
Third-party dependencies retain their own licenses. Publisher PDFs, third-party book illustrations and supplied original manuscripts are not redistributed.

> Uyulan, C. (2026). Thermal Preview Learning and Lyapunov Safeguarding for Satellite Appendage Vibration Control — reproducibility package (Version 2.0.0). Zenodo. https://doi.org/10.5281/zenodo.23126822

[CITATION.cff](CITATION.cff) provides machine-readable repository/archive citation metadata.
For scientific issues or proposed changes, follow [CONTRIBUTING.md](CONTRIBUTING.md).
