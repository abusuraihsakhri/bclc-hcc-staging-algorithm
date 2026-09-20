"""Stage-level treatment-pathway summaries for BCLC HCC staging."""
from typing import Any, Dict


def stratify_treatment(
    bclc_stage: str,
    milan_eligible: bool,
    child_pugh: str,
    ecog_ps: int,
    tumor_size_cm: float,
    tumor_count: int,
    portal_invasion: bool = False,
) -> Dict[str, Any]:
    """Return pathway-review flags without claiming definitive candidacy."""
    stage = str(bclc_stage).upper()
    if stage not in {"0", "A", "B", "C", "D"}:
        raise ValueError(f"Invalid BCLC stage: {bclc_stage!r}")

    return {
        "bclc_stage": stage,
        "treatment_pathway": _treatment_pathway(stage),
        "milan_criteria_eligible": bool(milan_eligible),
        "transplant_eligible": None,
        "transplant_evaluation_required": bool(milan_eligible or stage in {"0", "A", "B"}),
        "resection_review": stage in {"0", "A"},
        "locoregional_therapy_review": stage in {"0", "A", "B"},
        "systemic_therapy_review": stage in {"B", "C"},
        "supportive_care_review": stage == "D",
        "clinical_trial_review": stage in {"B", "C"},
        "input_context": {
            "child_pugh": child_pugh,
            "ecog_ps": ecog_ps,
            "tumor_size_cm": tumor_size_cm,
            "tumor_count": tumor_count,
            "portal_invasion": portal_invasion,
        },
        "limitations": (
            "Flags indicate topics for multidisciplinary review only. Definitive "
            "treatment or transplant eligibility requires additional variables."
        ),
    }


def _treatment_pathway(stage: str) -> str:
    return {
        "0": "Very-early-stage pathway: evaluate ablation or resection and liver/transplant context",
        "A": "Early-stage pathway: evaluate resection, ablation, and liver transplantation",
        "B": "Intermediate-stage pathway: subgroup for transplant, locoregional, or systemic treatment",
        "C": "Advanced-stage pathway: systemic therapy is the usual evidence-based first option",
        "D": "Supportive/palliative-care pathway; review transplant context when relevant",
    }[stage]
