# Barcelona Clinic Liver Cancer (BCLC) HCC Staging Algorithm

A clinically validated, pure Python staging and treatment allocation engine for hepatocellular carcinoma (HCC) based on the **2022 updated BCLC Staging and Treatment Strategy** (Reig et al., *Journal of Hepatology* 2022) and international EASL / AASLD clinical practice guidelines.

---

## Clinical Staging Architecture

The BCLC system integrates three critical prognostic pillars:
1. **Tumor Burden**: Number of nodules, maximum nodule diameter (cm), macrovascular portal vein invasion (PVI), and extrahepatic spread (EHS).
2. **Liver Function Reserve**: Child-Pugh class (A, B, or C).
3. **General Performance Status**: Eastern Cooperative Oncology Group (ECOG PS 0–4).

### Staging Rules & Treatment Allocations

| BCLC Stage | Category | Tumor Burden | Liver Function | ECOG PS | Primary Treatment Allocation | Median Survival | 5-Year Survival |
|:---|:---|:---|:---|:---|:---|:---|:---|
| **0** | **Very Early** | Solitary nodule $\le 2\text{ cm}$ | Child-Pugh A | 0 | Surgical Resection or Local Ablation (RFA/MWA) | $> 60\text{ months}$ | ~70% |
| **A** | **Early** | Solitary nodule (any size) OR $\le 3$ nodules each $\le 3\text{ cm}$ | Child-Pugh A–B | 0 | Liver Transplantation, Resection, or Ablation | $> 36\text{ months}$ | ~50% |
| **B** | **Intermediate** | Multifocal / multinodular ($>3$ nodules OR $>1$ nodule $>3\text{ cm}$) without vascular invasion | Child-Pugh A–B | 0 | Transarterial Chemoembolization (TACE) / systemic combination | ~20–30 months | ~25% |
| **C** | **Advanced** | Macrovascular invasion (PVI/vascular) OR Extrahepatic spread (N1/M1) | Child-Pugh A–B | 1–2 | First-line Systemic Therapy (Atezolizumab + Bevacizumab, Tremelimumab + Durvalumab, Sorafenib, Lenvatinib) | ~12–19 months | ~10% |
| **D** | **Terminal** | Any tumor burden | Child-Pugh C | 3–4 | Best Supportive Care (Palliative / Hospice); Transplant only if Milan-eligible with sole Child-Pugh C criterion | $< 3\text{ months}$ | ~0% |

---

## Milan Criteria for Liver Transplantation

The engine automatically evaluates eligibility for cadaveric or living-donor liver transplantation under the Mazzaferro (1996) **Milan Criteria**:
- **Single nodule**: Diameter $\le 5.0\text{ cm}$
- **Multiple nodules**: $\le 3$ nodules, each with diameter $\le 3.0\text{ cm}$
- **Exclusions**: Absence of macrovascular portal vein invasion and absence of extrahepatic metastases

---

## Features

- **2022 BCLC Precision:** Full support for very early (0), early (A), intermediate (B), advanced (C), and terminal (D) staging.
- **Milan Transplantation Calculator:** Instant organ allocation candidacy verification.
- **Batch CSV Processing:** High-throughput batch triage for multidisciplinary liver tumor boards.
- **Zero Runtime Dependencies:** Standalone implementation utilizing the Python Standard Library only.

---

## Installation & Requirements

- Python 3.10+ (tested on 3.10, 3.11, 3.12)
- Zero external runtime dependencies.

```bash
git clone https://github.com/abusuraihsakhri/bclc-hcc-staging-algorithm.git
cd bclc-hcc-staging-algorithm
```

---

## CLI Usage

### 1. Stage a Single HCC Patient
```bash
python cli.py single --tumor-count 1 --tumor-size-cm 1.8 --child-pugh-class A --ecog-ps 0
```

### 2. Evaluate an Advanced Case
```bash
python cli.py single --tumor-count 2 --tumor-size-cm 4.5 --child-pugh-class B --ecog-ps 1 --portal-vein-invasion
```

### 3. Batch Process Patient Cohorts from CSV
```bash
python cli.py batch --input sample.csv --output results.csv
```

---

## Python API Quickstart

```python
from bclc_staging import stage_bclc

# Stage patient with multifocal HCC
result = stage_bclc(
    tumor_count=3,
    tumor_size_cm=2.5,
    child_pugh_class="A",
    ecog_ps=0,
    portal_vein_invasion=False,
    extrahepatic_spread=False
)

print(f"BCLC Stage: {result['bclc_stage']} ({result['stage_name']})")
print(f"Allocation: {result['treatment_allocation']}")
print(f"Milan Eligible: {result['milan_criteria_eligible']}")
print(f"5-Year Survival: {result['five_year_survival_pct']}%")
```

---

## Testing & Verification

Run the comprehensive unit test suite:

```bash
python -m pytest -p no:zarr
```

