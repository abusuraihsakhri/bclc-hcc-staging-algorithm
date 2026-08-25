# BCLC (Barcelona Clinic Liver Cancer) Staging Algorithm

Real implementation of the BCLC staging system for hepatocellular carcinoma (HCC), per EASL-AASLD guidelines.

## What It Does

Stages HCC into five groups and provides treatment allocation:

| Stage | Name | Criteria | Treatment | 5-Year OS |
|-------|------|----------|-----------|-----------|
| 0 | Very Early | Single ≤2cm, Child-Pugh A, ECOG 0 | Ablation/Resection | ~70% |
| A | Early | Single or ≤3 nodules ≤3cm, CP A-B, ECOG 0 | Resection/Transplant/Ablation | ~50% |
| B | Intermediate | Multinodular, CP A-B, ECOG 0 | TACE | ~25% |
| C | Advanced | Portal invasion/N1/M1, CP A-B, ECOG 1-2 | Systemic therapy | ~10% |
| D | Terminal | Any tumor, Child-Pugh C or ECOG 3-4 | Best supportive care | ~0% |

Also evaluates **Milan criteria** eligibility and **transplant eligibility**.

## Installation

Zero dependencies — Python 3.7+ stdlib only.

## Usage

### Single Patient

```bash
python bclc_staging.py single \
  --tumor-count 1 \
  --tumor-size-cm 2.5 \
  --child-pugh-class A \
  --ecog-ps 0
```

With invasion flags:
```bash
python bclc_staging.py single \
  --tumor-count 3 --tumor-size-cm 3.0 \
  --child-pugh-class B --ecog-ps 0 \
  --portal-vein-invasion
```

### Batch Processing

```bash
python bclc_staging.py batch -i patients.csv -o results.csv
```

CSV columns: `tumor_count`, `tumor_size_cm`, `child_pugh_class`, `ecog_ps`, `portal_vein_invasion`, `extrahepatic_spread`, `lymph_node_metastasis`, `vascular_invasion`

### Python API

```python
from bclc_staging import stage_bclc

result = stage_bclc(
    tumor_count=1, tumor_size_cm=2.0,
    child_pugh_class="A", ecog_ps=0,
)
print(result["bclc_stage"])              # "0"
print(result["treatment_allocation"])     # "Ablation or Resection"
print(result["milan_criteria_eligible"])  # True
```

## Running Tests

```bash
python -m pytest test_bclc_staging.py -v
```

## Clinical Reference

Llovet JM et al. Hepatocellular carcinoma. Lancet. 2021;397(10287):1849-1862.

European Association for the Study of the Liver. EASL Clinical Practice Guidelines: Management of hepatocellular carcinoma. J Hepatol. 2018;69(1):182-236.

## License

MIT
