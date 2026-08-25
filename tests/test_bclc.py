"""Tests for bclc-hcc-staging-algorithm core and enrichment."""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


# --- Core Staging Tests ---
def test_very_early_stage():
    from bclc_staging import calculate_metrics
    r = calculate_metrics(tumor_size_cm=1.5, tumor_count=1, ecog_ps=0, child_pugh="A",
                          vascular_invasion=False, extrahepatic_spread=False)
    assert r["bclc_stage"] == "0"


def test_early_stage_milan():
    from bclc_staging import calculate_metrics
    r = calculate_metrics(tumor_size_cm=2.5, tumor_count=2, ecog_ps=0, child_pugh="A",
                          vascular_invasion=False, extrahepatic_spread=False)
    assert r["bclc_stage"] == "A"
    assert r["milan_criteria_eligible"] is True


def test_intermediate_stage():
    from bclc_staging import calculate_metrics
    r = calculate_metrics(tumor_size_cm=6.0, tumor_count=1, ecog_ps=0, child_pugh="A",
                          vascular_invasion=False, extrahepatic_spread=False)
    assert r["bclc_stage"] == "B"


def test_advanced_stage():
    from bclc_staging import calculate_metrics
    r = calculate_metrics(tumor_size_cm=3.0, tumor_count=1, ecog_ps=1, child_pugh="A",
                          portal_invasion=True, extrahepatic_spread=False)
    assert r["bclc_stage"] == "C"


def test_terminal_stage_child_c():
    from bclc_staging import calculate_metrics
    r = calculate_metrics(tumor_size_cm=3.0, tumor_count=1, ecog_ps=0, child_pugh="C")
    assert r["bclc_stage"] == "D"


def test_terminal_stage_ecog():
    from bclc_staging import calculate_metrics
    r = calculate_metrics(tumor_size_cm=3.0, tumor_count=1, ecog_ps=3, child_pugh="A")
    assert r["bclc_stage"] == "D"


def test_milan_outside():
    from bclc_staging import calculate_metrics
    r = calculate_metrics(tumor_size_cm=6.0, tumor_count=1, ecog_ps=0, child_pugh="A")
    assert r["milan_criteria_eligible"] is False


def test_intermediate_with_b_child():
    from bclc_staging import calculate_metrics
    r = calculate_metrics(tumor_size_cm=4.0, tumor_count=1, ecog_ps=0, child_pugh="B")
    assert r["bclc_stage"] == "B"


# --- Serial Tracking Tests ---
def test_bclc_tracker_timeline():
    from serial_tracking import BCLCTracker
    t = BCLCTracker()
    t.add_staging("P001", "2026-01-01", "A")
    t.add_staging("P001", "2026-06-01", "B")
    t.add_staging("P001", "2026-12-01", "C")
    tl = t.get_timeline("P001")
    assert len(tl) == 3
    assert [r.bclc_stage for r in tl] == ["A", "B", "C"]


def test_bclc_tracker_progression():
    from serial_tracking import BCLCTracker
    t = BCLCTracker()
    t.add_staging("P001", "2026-01-01", "A")
    t.add_staging("P001", "2026-06-01", "B")
    alerts = t.detect_progression("P001")
    assert len(alerts) == 1
    assert alerts[0]["type"] == "BCLC_PROGRESSION"


def test_bclc_tracker_downstaging():
    from serial_tracking import BCLCTracker
    t = BCLCTracker()
    t.add_staging("P001", "2026-01-01", "B")
    t.add_staging("P001", "2026-06-01", "A")
    alerts = t.detect_progression("P001")
    assert alerts[0]["type"] == "BCLC_DOWNSTAGING"


def test_bclc_tracker_cohort():
    from serial_tracking import BCLCTracker
    t = BCLCTracker()
    t.add_staging("P001", "2026-01-01", "A")
    t.add_staging("P002", "2026-01-01", "C")
    t.add_staging("P003", "2026-01-01", "A")
    s = t.get_cohort_summary()
    assert s["total"] == 3
    assert s["stage_distribution"]["A"] == 2


# --- Alert Escalation Tests ---
def test_ecog_decline_alert():
    from alert_escalation import BCLCAlertEngine
    e = BCLCAlertEngine()
    alerts = e.evaluate("P001", "A", True, ecog_ps=2, child_pugh="A", tumor_count=1, tumor_size_cm=2.0)
    assert any(a["alert_type"] == "ECOG_DECLINE" for a in alerts)


def test_transplant_eligible_alert():
    from alert_escalation import BCLCAlertEngine
    e = BCLCAlertEngine()
    alerts = e.evaluate("P001", "B", True, ecog_ps=0, child_pugh="A", tumor_count=2, tumor_size_cm=3.0)
    assert any(a["alert_type"] == "TRANSPLANT_ELIGIBLE" for a in alerts)


def test_terminal_alert():
    from alert_escalation import BCLCAlertEngine
    e = BCLCAlertEngine()
    alerts = e.evaluate("P001", "D", False, ecog_ps=3, child_pugh="C", tumor_count=1, tumor_size_cm=3.0)
    assert any(a["alert_type"] == "TERMINAL_STAGE" for a in alerts)


def test_child_c_alert():
    from alert_escalation import BCLCAlertEngine
    e = BCLCAlertEngine()
    alerts = e.evaluate("P001", "B", False, ecog_ps=0, child_pugh="C", tumor_count=1, tumor_size_cm=3.0)
    assert any(a["alert_type"] == "DECOMPENSATED_CIRRHOSIS" for a in alerts)


# --- Stratification Tests ---
def test_stratify_resectable():
    from stratification import stratify_treatment
    r = stratify_treatment("0", True, "A", 0, 1.5, 1)
    assert r["resectable"] is True
    assert r["best_supportive_care"] is False


def test_stratify_tace():
    from stratification import stratify_treatment
    r = stratify_treatment("B", True, "A", 0, 6.0, 4)
    assert r["tace_candidate"] is True
    assert r["transplant_eligible"] is True  # milan + good liver


def test_stratify_systemic():
    from stratification import stratify_treatment
    r = stratify_treatment("C", False, "A", 1, 3.0, 1, portal_invasion=True)
    assert r["systemic_therapy"] is True


def test_stratify_bsc():
    from stratification import stratify_treatment
    r = stratify_treatment("D", False, "C", 3, 3.0, 1)
    assert r["best_supportive_care"] is True


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            failed += 1
            print(f"FAIL: {t.__name__}: {e}")
    print(f"\n{passed} passed, {failed} failed out of {len(tests)}")
