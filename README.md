# Bclc Hcc Staging Algorithm

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Alert escalation engine for BCLC staging.

Barcelona Clinic Liver Cancer (BCLC) Staging Algorithm

Stages hepatocellular carcinoma (HCC) into stages 0, A, B, C, D based on
tumor characteristics, liver function (Child-Pugh), and performance status (ECOG).
Provides treatment allocation per EASL-AASLD guidelines.

Staging criteria:
  Stage 0 (Very Early): Single ≤2cm, Child-Pugh A, ECOG 0
  Stage A (Early):      Single or ≤3 nodules each ≤3cm, Child-Pugh A-B, ECOG 0
  Stage B (Intermediate): Multinodular, Child-Pugh A-B, ECOG 0
  Stage C (Advanced):   Portal invasion/N1/M1, Child-Pugh A-B, ECOG 1-2
  Stage D (Terminal):   Any tumor, Child-Pugh C, ECOG 3-4

Treatment by stage:
  0 → Ablation/Resection
  A → Resection / Transplant / Ablation
  B → TACE (Transarterial chemoembolization)
  C → Systemic therapy
  D → Best supportive care

Zero-dependency Python implementation.
License: MIT

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`BCLCAlertEngine`** — dedicated module for b c l c alert engine evaluation and state verification.
- **`BCLCStagingRecord`** — dedicated module for b c l c staging record evaluation and state verification.
- **`BCLCTracker`** — dedicated module for b c l c tracker evaluation and state verification.
- **`HCCAssessment`**: One timepoint of HCC tumor burden + liver function.
- **`StagingEpisode`** — dedicated module for staging episode evaluation and state verification.

---

## 📐 Mathematical Formulation & Logic

```text
  return (
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --input data.csv
```

### Parameter Reference
- `--interactive`: Launch guided terminal interactive wizard.
- `--input <path>`: Evaluate input from JSON or CSV specification.
- `--json`: Output deterministic structured results in JSON format.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `Patient_ID` | Parameter / observation metric | Required |
| `v1` | Parameter / observation metric | Required |
| `v2` | Parameter / observation metric | Required |
| `v3` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t bclc-hcc-staging-algorithm .
docker run -p 8000:8000 bclc-hcc-staging-algorithm
```
