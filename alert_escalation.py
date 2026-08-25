"""Alert escalation engine for BCLC staging."""
from typing import List


class BCLCAlertEngine:
    def __init__(self):
        self.active_alerts: List[dict] = []

    def evaluate(self, patient_id: str, bclc_stage: str, milan_eligible: bool,
                 ecog_ps: int, child_pugh: str, tumor_count: int,
                 tumor_size_cm: float, **kwargs) -> List[dict]:
        alerts = []
        if ecog_ps >= 2:
            alerts.append({
                "alert_type": "ECOG_DECLINE",
                "severity": "CRITICAL",
                "patient_id": patient_id,
                "message": f"ECOG PS {ecog_ps}: patient no longer eligible for curative therapy",
                "routed_to": "Palliative care / tumor board",
            })
        if bclc_stage == "D":
            alerts.append({
                "alert_type": "TERMINAL_STAGE",
                "severity": "CRITICAL",
                "patient_id": patient_id,
                "message": "BCLC-D: best supportive care indicated",
                "routed_to": "Palliative care / hospice",
            })
        if bclc_stage == "B" and milan_eligible:
            alerts.append({
                "alert_type": "TRANSPLANT_ELIGIBLE",
                "severity": "HIGH",
                "patient_id": patient_id,
                "message": "BCLC-B with Milan criteria: consider TACE downstaging to transplant",
                "routed_to": "Transplant hepatology",
            })
        if bclc_stage == "C":
            alerts.append({
                "alert_type": "ADVANCED_HCC",
                "severity": "HIGH",
                "patient_id": patient_id,
                "message": "BCLC-C: systemic therapy required, consider clinical trials",
                "routed_to": "Oncology",
            })
        if tumor_count > 3 or tumor_size_cm > 5.0:
            alerts.append({
                "alert_type": "ADVANCED_TUMOR_BURDEN",
                "severity": "MEDIUM",
                "patient_id": patient_id,
                "message": f"Tumor burden: {tumor_count} nodules, {tumor_size_cm}cm. Re-evaluate staging.",
                "routed_to": "Hepatology / tumor board",
            })
        if child_pugh == "C":
            alerts.append({
                "alert_type": "DECOMPENSATED_CIRRHOSIS",
                "severity": "CRITICAL",
                "patient_id": patient_id,
                "message": "Child-Pugh C: liver function precludes curative therapy",
                "routed_to": "Transplant evaluation / palliative care",
            })
        self.active_alerts.extend(alerts)
        return alerts

    def get_active_alerts(self, patient_id: str = None) -> List[dict]:
        if patient_id:
            return [a for a in self.active_alerts if a["patient_id"] == patient_id]
        return self.active_alerts

    def dismiss_alert(self, alert_type: str, patient_id: str) -> bool:
        for i, a in enumerate(self.active_alerts):
            if a["alert_type"] == alert_type and a["patient_id"] == patient_id:
                self.active_alerts.pop(i)
                return True
        return False
