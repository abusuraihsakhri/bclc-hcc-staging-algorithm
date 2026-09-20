# BCLC HCC Staging Algorithm

### [Open the Live Application →](https://abusuraihsakhri.github.io/bclc-hcc-staging-algorithm/)

A zero-runtime-dependency Python implementation of the **Barcelona Clinic Liver Cancer (BCLC) 2026** stage-defining framework for confirmed hepatocellular carcinoma (HCC), with CLI, CSV batch processing, longitudinal helpers, regression tests, and a browser interface powered by Pyodide.

> **Scope:** education, research, reproducible data processing, and software validation. This repository does not replace multidisciplinary clinical assessment, transplant-center criteria, or current treatment guidelines.

## What it implements

The core engine uses tumor burden, macrovascular invasion / metastatic spread, HCC-attributable ECOG performance status, and decompensation/transplant context to assign BCLC stage 0, A, B, C, or D.

Important 2026 boundaries implemented in the test suite include:

- **BCLC 0:** one tumor up to 2 cm, without advanced-stage features.
- **BCLC A:** one tumor over 2 cm with no upper size cutoff, or up to 3 tumors with none over 3 cm.
- **BCLC B:** multifocal disease with more than 3 tumors, or 2–3 tumors with at least one over 3 cm, without BCLC-C/D features.
- **BCLC C:** macrovascular invasion, extrahepatic/metastatic spread, or HCC-attributable ECOG 1–2.
- **BCLC D:** severe HCC-attributable performance-status impairment or decompensated liver disease when transplantation is not an option.

Child-Pugh class is retained as descriptive/compatibility input, but it is not used alone to force a BCLC stage.

## Features

- BCLC 2026 staging with explicit boundary validation.
- Milan tumor-burden screening with macrovascular/metastatic exclusions.
- Transplant candidacy kept separate from Milan criteria.
- Single-case JSON CLI output.
- CSV batch processing with row-level validation errors.
- Longitudinal stage tracking and current OPTN HCC downstaging screening helper.
- Dependency-free stress/smoke simulator.
- Modern light-first browser UI with an explicit persistent dark-mode option.
- Static browser UI that executes the same Python module through Pyodide.
- GitHub Actions regression matrix for Python 3.10–3.14 and a headless-browser Pyodide smoke test.

## CLI

No runtime package installation is required.

```bash
python cli.py single \
  --tumor-count 1 \
  --tumor-size-cm 6 \
  --child-pugh-class A \
  --ecog-ps 0
```

Batch processing:

```bash
python cli.py batch --input sample.csv --output results.csv
```

Use `python cli.py single --help` for decompensation, transplant-context, invasion, spread, and ECOG-attribution flags.

## Python API

```python
from bclc_staging import stage_bclc

result = stage_bclc(
    tumor_count=3,
    tumor_size_cm=2.5,
    child_pugh_class="A",
    ecog_ps=0,
    portal_vein_invasion=False,
    extrahepatic_spread=False,
)

print(result["classification"])
print(result["milan_criteria_eligible"])
```

## Browser application

`index.html` provides a responsive light-first interface with an explicit light/dark theme toggle, compact stage-defining inputs, accessible focus states, and a sticky result panel on larger screens. It loads Pyodide 0.29.5 and imports the repository's `bclc_staging.py` directly, so the browser UI does not maintain a separate JavaScript copy of the staging rules.

Case inputs are not stored or submitted by the application. The selected UI theme may be stored locally in the browser. Initial page load downloads the Pyodide runtime from jsDelivr; calculations then execute in the browser.

## Testing

```bash
python -m pip install pytest
python -m compileall -q .
python -m pytest -q
python simulator.py 25
```

CI additionally exercises the CLI, sample CSV workflow, Python 3.10–3.14, and the browser/Pyodide path in headless Chrome.

## Clinical and technical limitations

- Assumes HCC has already been established; it is not a diagnostic model.
- ECOG-based upstaging should reflect symptoms attributable to HCC rather than unrelated disability.
- `vascular_invasion` is treated as **macrovascular** invasion for backward compatibility.
- Milan criteria do not establish transplant eligibility by themselves.
- OPTN downstaging output is a screening aid; current OPTN policy and transplant-center review are authoritative.
- Treatment strings are stage-level pathway summaries, not patient-specific prescriptions.
- No outcome model or individualized survival prediction is implemented.

## References

1. Reig M, et al. BCLC strategy for prognosis prediction and treatment recommendations: The 2026 update. *Journal of Hepatology*. 2026;84(3):631-654. doi:10.1016/j.jhep.2025.10.020
2. OPTN Policies, Policy 9.5.I: Requirements for Hepatocellular Carcinoma (HCC) MELD or PELD Score Exceptions. Current policy should be checked before clinical or allocation use.

## License

MIT License. See [LICENSE](LICENSE).
