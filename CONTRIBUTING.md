# Contributing

Report numerical or documentation issues through this repository's Issues page. Include the script, software versions,
case identifier or random seed, exact reproduction command, observed result and expected result. Attach a small numerical
example when possible. Distinguish implementation errors from the declared model and validation limitations.

For a pull request, explain the scientific effect of the change and compare the relevant baseline before updating a report
or figure. Preserve paired case seeds and common metric definitions. Run `python code/verify_package.py` before modifying
the frozen payload. Regenerated data must be identified as changed results and receive an updated manifest; do not silently
replace published benchmarks. Run numerical checks appropriate to the changed equations or controller.

Source code contributions use MIT. Generated data, figures and documentation use CC BY 4.0, as documented in README and
the license files. Preserve third-party attribution. Do not add licensed publisher PDFs or book illustrations.

Scientific contact: Caglar Uyulan, caglar.uyulan@ikcu.edu.tr.
