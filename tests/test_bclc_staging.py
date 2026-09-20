import csv
import pytest

from bclc_staging import (
    _check_milan,
    _parse_bool,
    calculate_metrics,
    process_batch,
    stage_bclc,
)


@pytest.mark.parametrize(
    ("count", "size", "stage"),
    [
        (1, 1.5, "0"),
        (1, 2.0, "0"),
        (1, 2.01, "A"),
        (1, 6.0, "A"),
        (2, 3.0, "A"),
        (3, 3.0, "A"),
        (2, 3.1, "B"),
        (4, 2.0, "B"),
    ],
)
def test_tumor_burden_boundaries(count, size, stage):
    result = stage_bclc(tumor_count=count, tumor_size_cm=size)
    assert result["bclc_stage"] == stage


@pytest.mark.parametrize(
    "kwargs",
    [
        {"portal_vein_invasion": True},
        {"vascular_invasion": True},
        {"extrahepatic_spread": True},
        {"lymph_node_metastasis": True},
    ],
)
def test_advanced_features_are_stage_c(kwargs):
    result = stage_bclc(tumor_count=1, tumor_size_cm=1.5, **kwargs)
    assert result["bclc_stage"] == "C"


@pytest.mark.parametrize("ecog", [1, 2])
def test_hcc_related_ecog_1_2_is_stage_c(ecog):
    assert stage_bclc(tumor_count=1, tumor_size_cm=3, ecog_ps=ecog)["bclc_stage"] == "C"


def test_hcc_related_ecog_3_is_stage_d():
    assert stage_bclc(tumor_count=1, tumor_size_cm=3, ecog_ps=3)["bclc_stage"] == "D"


def test_unrelated_ecog_does_not_upstage():
    result = stage_bclc(
        tumor_count=1,
        tumor_size_cm=3,
        ecog_ps=3,
        ecog_cancer_related=False,
    )
    assert result["bclc_stage"] == "A"


def test_child_pugh_c_alone_does_not_force_stage_d():
    result = stage_bclc(
        tumor_count=1,
        tumor_size_cm=3,
        child_pugh_class="C",
        ecog_ps=0,
    )
    assert result["bclc_stage"] == "A"


def test_decompensation_requires_transplant_context():
    with pytest.raises(ValueError, match="transplant_candidate"):
        stage_bclc(
            tumor_count=1,
            tumor_size_cm=3,
            liver_decompensation=True,
        )


def test_decompensated_not_transplant_candidate_is_stage_d():
    result = stage_bclc(
        tumor_count=1,
        tumor_size_cm=3,
        liver_decompensation=True,
        transplant_candidate=False,
    )
    assert result["bclc_stage"] == "D"


def test_decompensated_transplant_candidate_not_automatically_d():
    result = stage_bclc(
        tumor_count=1,
        tumor_size_cm=3,
        liver_decompensation=True,
        transplant_candidate=True,
    )
    assert result["bclc_stage"] == "A"


def test_milan_tumor_burden_and_exclusions():
    assert _check_milan(1, 5.0)
    assert _check_milan(3, 3.0)
    assert not _check_milan(1, 5.1)
    assert not _check_milan(4, 2.0)
    assert not _check_milan(1, 3.0, macrovascular_invasion=True)
    assert not _check_milan(1, 3.0, extrahepatic_spread=True)
    assert not _check_milan(1, 3.0, lymph_node_metastasis=True)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"tumor_count": 0, "tumor_size_cm": 2},
        {"tumor_count": 1, "tumor_size_cm": 0},
        {"tumor_count": 1, "tumor_size_cm": -1},
        {"tumor_count": 1, "tumor_size_cm": 2, "ecog_ps": 5},
        {"tumor_count": 1, "tumor_size_cm": 2, "child_pugh_class": "D"},
    ],
)
def test_invalid_inputs(kwargs):
    with pytest.raises(ValueError):
        stage_bclc(**kwargs)


def test_boolean_parser_is_strict():
    assert _parse_bool("YES") is True
    assert _parse_bool("0") is False
    with pytest.raises(ValueError):
        _parse_bool("maybe")


def test_compatibility_wrapper():
    result = calculate_metrics(
        tumor_count=1,
        tumor_size_cm=6.0,
        child_pugh="A",
        portal_invasion=False,
    )
    assert result["bclc_stage"] == "A"


def test_result_does_not_invent_transplant_eligibility():
    result = stage_bclc(tumor_count=1, tumor_size_cm=3)
    assert result["milan_criteria_eligible"] is True
    assert result["transplant_eligible"] is None
    assert result["five_year_survival_pct"] is None


def test_batch_reports_row_level_validation_errors(tmp_path):
    source = tmp_path / "input.csv"
    output = tmp_path / "output.csv"
    source.write_text(
        "tumor_count,tumor_size_cm,child_pugh_class,ecog_ps,portal_vein_invasion\n"
        "1,6,A,0,false\n"
        "2,3,A,0,maybe\n",
        encoding="utf-8",
    )
    assert process_batch(str(source), str(output)) == 2
    with output.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["bclc_stage"] == "A"
    assert rows[0]["error"] == ""
    assert rows[1]["bclc_stage"] == ""
    assert "Invalid boolean" in rows[1]["error"]
