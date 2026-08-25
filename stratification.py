"""Patient stratification for BCLC HCC staging."""
from typing import Dict, Any


def stratify_treatment(bclc_stage: str, milan_eligible: bool, child_pugh: str,
                       ecog_ps: int, tumor_size_cm: float, tumor_count: int,
                       portal_invasion: bool = False) -> Dict[str, Any]:
    """Stratify patient into BCLC treatment pathway."""
    treatment_pathway = _treatment_pathway(bclc_stage, milan_eligible, child_pugh)
    clinical_trials = _trial_eligibility(bclc_stage, ecog_ps, child_pugh)
    eligibility_checklist = _eligibility_checklist(bclc_stage, milan_eligible, child_pugh, ecog_ps)

    return {
        "bclc_stage": bclc_stage,
        "treatment_pathway": treatment_pathway,
        "transplant_eligible": milan_eligible and child_pugh in ("A", "B") and ecog_ps <= 1,
        "resectable": bclc_stage in ("0", "A") and child_pugh == "A" and ecog_ps == 0,
        "tace_candidate": bclc_stage == "B" and child_pugh in ("A", "B") and ecog_ps <= 2,
        "systemic_therapy": bclc_stage == "C",
        "best_supportive_care": bclc_stage == "D",
        "clinical_trial_eligible": clinical_trials,
        "eligibility_checklist": eligibility_checklist,
    }


def _treatment_pathway(stage, milan, child):
    if stage == "D":
        return "Best supportive care, palliative symptom management, hospice referral"
    if stage == "C":
        return "First-line: atezolizumab + bevacizumab. Second-line: sorafenib/lenvatinib/ramucirumab"
    if stage == "B":
        if milan:
            return "TACE → reassess for transplant if downstaged to within Milan"
        return "TACE (standard of care). Repeat TACE q6-8 weeks. Assess for transplant eligibility."
    if stage == "A":
        if milan:
            return "Liver transplant (preferred) or surgical resection if adequate liver reserve"
        return "Surgical resection. Consider living donor transplant if available."
    return "Curative resection or local ablation (RFA/MWA)"


def _trial_eligibility(stage, ecog, child):
    if child == "C" or ecog >= 3:
        return ["Palliative-only trials"]
    trials = []
    if stage in ("C", "B"):
        trials.append("Immunotherapy combination trials")
        trials.append("Targeted therapy trials")
    if stage in ("0", "A", "B"):
        trials.append("Adjuvant therapy trials")
        trials.append("Locoregional therapy trials")
    return trials


def _eligibility_checklist(stage, milan, child, ecog):
    checklist = []
    if stage in ("0", "A"):
        checklist.append(("CT/MRI staging", True))
        checklist.append(("Liver function assessment", child == "A"))
        checklist.append(("ECOG performance status", ecog <= 1))
        checklist.append(("Milan criteria", milan))
        checklist.append(("Vascular invasion assessment", True))
    elif stage == "B":
        checklist.append(("TACE candidacy", child in ("A", "B")))
        checklist.append(("Transplant evaluation", milan))
        checklist.append(("ECOG assessment", ecog <= 2))
    elif stage == "C":
        checklist.append(("Systemic therapy eligibility", True))
        checklist.append(("Clinical trial screening", True))
        checklist.append(("Liver function preserved", child in ("A", "B")))
    else:
        checklist.append(("Palliative care referral", True))
        checklist.append(("Symptom management plan", True))
    return checklist
