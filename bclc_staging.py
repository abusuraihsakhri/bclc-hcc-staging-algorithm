#!/usr/bin/env python3
"""
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
"""

import argparse
import csv
import json
import sys
from typing import Dict, Any, Optional


# ---------------------------------------------------------------------------
# BCLC stage definitions
# ---------------------------------------------------------------------------

BCLC_STAGES = {
    "0": {
        "name": "Very Early",
        "treatment": "Ablation or Resection",
        "median_survival_months": ">60",
        "five_year_survival_pct": 70.0,
    },
    "A": {
        "name": "Early",
        "treatment": "Resection / Liver Transplant / Ablation",
        "median_survival_months": ">36",
        "five_year_survival_pct": 50.0,
    },
    "B": {
        "name": "Intermediate",
        "treatment": "TACE (Transarterial Chemoembolization)",
        "median_survival_months": "20",
        "five_year_survival_pct": 25.0,
    },
    "C": {
        "name": "Advanced",
        "treatment": "Systemic Therapy (Atezolizumab+Bevacizumab, Sorafenib, Lenvatinib)",
        "median_survival_months": "12",
        "five_year_survival_pct": 10.0,
    },
    "D": {
        "name": "Terminal",
        "treatment": "Best Supportive Care",
        "median_survival_months": "<3",
        "five_year_survival_pct": 0.0,
    },
}

# Milan criteria for transplant eligibility
MILAN_SINGLE_MAX_CM = 5.0
MILAN_MULTI_MAX_COUNT = 3
MILAN_MULTI_MAX_CM = 3.0


# ---------------------------------------------------------------------------
# Core staging logic
# ---------------------------------------------------------------------------

def stage_bclc(
    tumor_count: int = 1,
    tumor_size_cm: float = 0.0,
    child_pugh_class: str = "A",
    ecog_ps: int = 0,
    portal_vein_invasion: bool = False,
    extrahepatic_spread: bool = False,
    lymph_node_metastasis: bool = False,
    vascular_invasion: bool = False,
) -> Dict[str, Any]:
    """
    Stage HCC using the BCLC algorithm.

    Parameters:
        tumor_count: Number of tumor nodules (default 1)
        tumor_size_cm: Largest tumor diameter in cm
        child_pugh_class: 'A', 'B', or 'C'
        ecog_ps: ECOG performance status 0-4
        portal_vein_invasion: Portal vein invasion present
        extrahepatic_spread: Extrahepatic metastases (M1)
        lymph_node_metastasis: Regional/distant lymph node involvement (N1)
        vascular_invasion: Micro or macro vascular invasion

    Returns:
        Dict with BCLC stage, treatment, survival estimates, and details.
    """
    # Validate inputs
    cp = child_pugh_class.strip().upper()
    if cp not in ("A", "B", "C"):
        raise ValueError(f"Invalid Child-Pugh class '{child_pugh_class}'. Use A, B, or C.")
    if ecog_ps < 0 or ecog_ps > 4:
        raise ValueError(f"ECOG PS must be 0-4, got {ecog_ps}")
    if tumor_count < 0:
        raise ValueError("Tumor count must be non-negative")
    if tumor_size_cm < 0:
        raise ValueError("Tumor size must be non-negative")

    # Determine stage
    stage = _determine_stage(
        tumor_count, tumor_size_cm, cp, ecog_ps,
        portal_vein_invasion, extrahepatic_spread,
        lymph_node_metastasis, vascular_invasion,
    )

    stage_info = BCLC_STAGES[stage]

    # Milan criteria
    milan_eligible = _check_milan(tumor_count, tumor_size_cm)

    # Transplant eligibility (Milan + adequate liver function)
    transplant_eligible = milan_eligible and cp in ("A", "B") and ecog_ps <= 1

    return {
        "tool": "bclc-hcc-staging-algorithm",
        "bclc_stage": stage,
        "stage_name": stage_info["name"],
        "treatment_allocation": stage_info["treatment"],
        "median_survival_months": stage_info["median_survival_months"],
        "five_year_survival_pct": stage_info["five_year_survival_pct"],
        "milan_criteria_eligible": milan_eligible,
        "transplant_eligible": transplant_eligible,
        "classification": f"BCLC Stage {stage} ({stage_info['name']})",
        "clinical_recommendation": _recommendation(stage, milan_eligible, cp, ecog_ps),
        "inputs": {
            "tumor_count": tumor_count,
            "tumor_size_cm": tumor_size_cm,
            "child_pugh_class": cp,
            "ecog_ps": ecog_ps,
            "portal_vein_invasion": portal_vein_invasion,
            "extrahepatic_spread": extrahepatic_spread,
            "lymph_node_metastasis": lymph_node_metastasis,
            "vascular_invasion": vascular_invasion,
        },
    }


def calculate_metrics(tumor_size_cm: float = 0.0, tumor_count: int = 1,
                      ecog_ps: int = 0, child_pugh: str = "A",
                      child_pugh_class: Optional[str] = None,
                      portal_invasion: bool = False,
                      portal_vein_invasion: Optional[bool] = None,
                      extrahepatic_spread: bool = False,
                      lymph_node_metastasis: bool = False,
                      vascular_invasion: bool = False, **kwargs) -> Dict[str, Any]:
    """Compatibility alias wrapping stage_bclc."""
    cp = child_pugh_class if child_pugh_class is not None else child_pugh
    pvi = portal_vein_invasion if portal_vein_invasion is not None else portal_invasion
    return stage_bclc(
        tumor_count=tumor_count,
        tumor_size_cm=tumor_size_cm,
        child_pugh_class=cp,
        ecog_ps=ecog_ps,
        portal_vein_invasion=pvi,
        extrahepatic_spread=extrahepatic_spread,
        lymph_node_metastasis=lymph_node_metastasis,
        vascular_invasion=vascular_invasion,
    )


def _determine_stage(
    count, size, cp, ecog, portal, extrahepatic, nodes, vascular
) -> str:
    """
    Apply the BCLC decision tree.

    Priority order (top-down):
    1. Child-Pugh C or ECOG 3-4 → Stage D (terminal)
    2. Macrovascular invasion, extrahepatic spread, or N1 → Stage C (advanced)
       (when ECOG 1-2 and Child-Pugh A-B)
    3. ECOG 1-2 with Child-Pugh A-B → Stage C (advanced)
    4. Multinodular or Child-Pugh B → Stage B (intermediate)
    5. Single ≤2cm, Child-Pugh A, ECOG 0 → Stage 0 (very early)
    6. Single or ≤3 nodules ≤3cm, Child-Pugh A-B, ECOG 0 → Stage A (early)
    7. Default → Stage B
    """
    # Stage D: terminal
    if cp == "C":
        return "D"
    if ecog >= 3:
        return "D"

    # Stage C: advanced - portal/hepatic vein invasion, extrahepatic, nodes
    has_major_invasion = portal or vascular
    has_metastasis = extrahepatic or nodes

    if has_major_invasion or has_metastasis:
        if cp in ("A", "B") and ecog <= 2:
            return "C"

    # ECOG 1-2 with preserved liver function → Stage C
    if ecog >= 1 and cp in ("A", "B"):
        return "C"

    # Stage B: intermediate - multinodular, or larger tumors, or Child-Pugh B
    if cp == "B":
        return "B"
    if count > 3:
        return "B"
    if count > 1 and size > 3.0:
        return "B"

    # Stage 0: very early - single ≤2cm, Child-Pugh A, ECOG 0
    if count == 1 and size <= 2.0 and cp == "A" and ecog == 0:
        return "0"

    # Stage A: early - single tumor or up to 3 nodules each ≤3cm
    if cp in ("A", "B") and ecog == 0:
        if count == 1 and size <= 5.0:
            return "A"
        if count <= 3 and size <= 3.0:
            return "A"

    # Default fallback
    return "B"


def _check_milan(count: int, size: float) -> bool:
    """Check if tumor meets Milan criteria for transplant."""
    if count == 1 and size <= MILAN_SINGLE_MAX_CM:
        return True
    if count <= MILAN_MULTI_MAX_COUNT and size <= MILAN_MULTI_MAX_CM:
        return True
    return False


def _recommendation(stage: str, milan: bool, cp: str, ecog: int) -> str:
    """Generate treatment recommendation."""
    if stage == "D":
        return (
            "Terminal stage. Best supportive care recommended. "
            "Hospice referral and symptom management. "
            "Transplant evaluation if Child-Pugh C is the sole criterion."
        )
    if stage == "C":
        return (
            "Advanced HCC. First-line systemic therapy: atezolizumab + bevacizumab "
            "(IMbrave150 regimen). Alternatives: sorafenib or lenvatinib. "
            "Consider clinical trial enrollment."
        )
    if stage == "B":
        if milan:
            return (
                "Intermediate HCC meeting Milan criteria. TACE as first-line therapy. "
                "Consider downstaging to transplant eligibility. "
                "Multidisciplinary tumor board evaluation recommended."
            )
        return (
            "Intermediate HCC outside Milan criteria. TACE recommended. "
            "Re-assess for transplant eligibility after downstaging. "
            "Consider combination therapy (TACE + systemic) for extensive disease."
        )
    if stage == "A":
        if milan:
            return (
                "Early HCC meeting Milan criteria. Liver transplantation is preferred "
                "(best long-term outcome). If not transplant candidate: surgical resection "
                "or local ablation (RFA/MWA)."
            )
        return (
            "Early HCC outside Milan criteria. Surgical resection if adequate liver "
            "reserve and no portal hypertension. Consider living donor transplant. "
            "Local ablation as alternative."
        )
    # Stage 0
    return (
        "Very early HCC. Surgical resection or local ablation (RFA/MWA) with "
        "curative intent. Excellent prognosis with 5-year survival ~70%. "
        "Standard post-treatment surveillance."
    )


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------

def process_batch(input_csv: str, output_csv: str) -> int:
    """Process a CSV of patients and write BCLC staging results."""
    with open(input_csv, mode="r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    out_fields = fieldnames + [
        "bclc_stage", "stage_name", "treatment_allocation",
        "milan_criteria_eligible", "transplant_eligible",
        "five_year_survival_pct", "clinical_recommendation",
    ]
    out_rows = []
    for r in rows:
        try:
            res = stage_bclc(
                tumor_count=int(r.get("tumor_count", 1)),
                tumor_size_cm=float(r.get("tumor_size_cm", 0)),
                child_pugh_class=r.get("child_pugh_class", "A"),
                ecog_ps=int(r.get("ecog_ps", 0)),
                portal_vein_invasion=_parse_bool(r.get("portal_vein_invasion", "false")),
                extrahepatic_spread=_parse_bool(r.get("extrahepatic_spread", "false")),
                lymph_node_metastasis=_parse_bool(r.get("lymph_node_metastasis", "false")),
                vascular_invasion=_parse_bool(r.get("vascular_invasion", "false")),
            )
            row_dict = dict(r)
            row_dict["bclc_stage"] = res["bclc_stage"]
            row_dict["stage_name"] = res["stage_name"]
            row_dict["treatment_allocation"] = res["treatment_allocation"]
            row_dict["milan_criteria_eligible"] = res["milan_criteria_eligible"]
            row_dict["transplant_eligible"] = res["transplant_eligible"]
            row_dict["five_year_survival_pct"] = res["five_year_survival_pct"]
            row_dict["clinical_recommendation"] = res["clinical_recommendation"]
        except (ValueError, KeyError) as e:
            row_dict = dict(r)
            row_dict["bclc_stage"] = f"ERROR: {e}"
            row_dict["stage_name"] = ""
            row_dict["treatment_allocation"] = ""
            row_dict["milan_criteria_eligible"] = ""
            row_dict["transplant_eligible"] = ""
            row_dict["five_year_survival_pct"] = ""
            row_dict["clinical_recommendation"] = ""
        out_rows.append(row_dict)

    with open(output_csv, mode="w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Processed {len(out_rows)} records -> {output_csv}")
    return len(out_rows)


def _parse_bool(val) -> bool:
    """Parse various boolean representations."""
    if isinstance(val, bool):
        return val
    s = str(val).strip().lower()
    return s in ("true", "1", "yes", "y", "t")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="BCLC (Barcelona Clinic Liver Cancer) Staging Algorithm"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Single evaluation
    sp = subparsers.add_parser("single", help="Stage a single HCC patient")
    sp.add_argument("--tumor-count", type=int, default=1,
                    help="Number of tumor nodules (default: 1)")
    sp.add_argument("--tumor-size-cm", type=float, default=0.0,
                    help="Largest tumor diameter in cm")
    sp.add_argument("--child-pugh-class", default="A",
                    choices=["A", "B", "C"],
                    help="Child-Pugh class (default: A)")
    sp.add_argument("--ecog-ps", type=int, default=0,
                    help="ECOG performance status 0-4 (default: 0)")
    sp.add_argument("--portal-vein-invasion", action="store_true",
                    help="Portal vein invasion present")
    sp.add_argument("--extrahepatic-spread", action="store_true",
                    help="Extrahepatic metastases present")
    sp.add_argument("--lymph-node-metastasis", action="store_true",
                    help="Lymph node metastasis present")
    sp.add_argument("--vascular-invasion", action="store_true",
                    help="Vascular invasion present")

    # Batch processing
    bp = subparsers.add_parser("batch", help="Batch process CSV file")
    bp.add_argument("-i", "--input", required=True, help="Input CSV file")
    bp.add_argument("-o", "--output", default="results.csv", help="Output CSV file")

    args = parser.parse_args(argv)

    if args.command == "single":
        result = stage_bclc(
            tumor_count=args.tumor_count,
            tumor_size_cm=args.tumor_size_cm,
            child_pugh_class=args.child_pugh_class,
            ecog_ps=args.ecog_ps,
            portal_vein_invasion=args.portal_vein_invasion,
            extrahepatic_spread=args.extrahepatic_spread,
            lymph_node_metastasis=args.lymph_node_metastasis,
            vascular_invasion=args.vascular_invasion,
        )
        print(json.dumps(result, indent=2))
    elif args.command == "batch":
        process_batch(args.input, args.output)


if __name__ == "__main__":
    main()
