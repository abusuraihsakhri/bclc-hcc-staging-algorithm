"""Alert helpers for longitudinal BCLC review."""
from typing import List


class BCLCAlertEngine:
    def __init__(self):
        self.active_alerts: List[dict] = []

    def evaluate(
        self,
        patient_id: str,
        bclc_stage: str,
        milan_eligible: bool,
        ecog_ps: int,
        child_pugh: str,
        tumor_count: int,
        tumor_size_cm: float,
        liver_decompensation: bool = False,
        transplant_candidate=None,
        **kwargs,
    ) -> List[dict]:
        stage = str(bclc_stage).upper()
        alerts = []

        if ecog_ps >= 3:
            alerts.append({
                "alert_type": "ECOG_SEVERE_IMPAIRMENT",
                "severity": "CRITICAL",
                "patient_id": patient_id,
                "message": f"ECOG PS {ecog_ps}: verify whether impairment is HCC-attributable and reassess BCLC stage.",
                "routed_to": "Multidisciplinary tumor board / supportive care",
            })
        elif ecog_ps in {1, 2}:
            alerts.append({
                "alert_type": "ECOG_DECLINE",
                "severity": "HIGH",
                "patient_id": patient_id,
                "message": f"ECOG PS {ecog_ps}: if HCC-attributable this is an advanced-stage feature.",
                "routed_to": "Multidisciplinary tumor board",
            })

        if stage == "D":
            alerts.append({
                "alert_type": "END_STAGE",
                "severity": "CRITICAL",
                "patient_id": patient_id,
                "message": "BCLC-D recorded: reassess reversible causes, transplant context, and supportive-care needs.",
                "routed_to": "Hepatology / palliative care",
            })
        elif stage == "C":
            alerts.append({
                "alert_type": "ADVANCED_HCC",
                "severity": "HIGH",
                "patient_id": patient_id,
                "message": "BCLC-C recorded: treatment plan requires advanced-stage multidisciplinary review.",
                "routed_to": "Hepatology / oncology",
            })

        if liver_decompensation and transplant_candidate is None:
            alerts.append({
                "alert_type": "TRANSPLANT_STATUS_REQUIRED",
                "severity": "CRITICAL",
                "patient_id": patient_id,
                "message": "Decompensation is present but transplant candidacy is unknown; BCLC-D cannot be inferred safely.",
                "routed_to": "Transplant hepatology",
            })

        if str(child_pugh).upper() == "C":
            alerts.append({
                "alert_type": "SEVERE_LIVER_DYSFUNCTION",
                "severity": "HIGH",
                "patient_id": patient_id,
                "message": "Child-Pugh C requires liver-function and transplant assessment; class alone does not define BCLC-D.",
                "routed_to": "Hepatology / transplant review",
            })

        if milan_eligible:
            alerts.append({
                "alert_type": "MILAN_TUMOR_BURDEN",
                "severity": "INFO",
                "patient_id": patient_id,
                "message": "Tumor burden is within Milan criteria; this does not by itself establish transplant eligibility.",
                "routed_to": "Transplant review when clinically appropriate",
            })

        if tumor_count > 3 or tumor_size_cm > 5.0:
            alerts.append({
                "alert_type": "HIGH_TUMOR_BURDEN",
                "severity": "MEDIUM",
                "patient_id": patient_id,
                "message": f"Tumor burden recorded as {tumor_count} nodule(s), largest {tumor_size_cm:g} cm.",
                "routed_to": "Multidisciplinary tumor board",
            })

        self.active_alerts.extend(alerts)
        return alerts

    def get_active_alerts(self, patient_id: str = None) -> List[dict]:
        if patient_id is None:
            return list(self.active_alerts)
        return [a for a in self.active_alerts if a["patient_id"] == patient_id]

    def dismiss_alert(self, alert_type: str, patient_id: str) -> bool:
        for index, alert in enumerate(self.active_alerts):
            if alert["alert_type"] == alert_type and alert["patient_id"] == patient_id:
                self.active_alerts.pop(index)
                return True
        return False
