import csv

from thznerve.runlog import RunLog


def test_runlog_appends_row(tmp_path):
    p = tmp_path / "log.csv"
    with RunLog("simX", "desc", path=p) as log:
        log.dof = 123
        log.outcome = "ok"
    with RunLog("simY", "desc2", path=p) as log:
        pass
    rows = list(csv.DictReader(p.open(encoding="utf-8")))
    assert [r["sim"] for r in rows] == ["simX", "simY"]
    assert rows[0]["dof"] == "123"
    assert float(rows[0]["peak_mem_gb"]) > 0


def test_runlog_records_failure(tmp_path):
    p = tmp_path / "log.csv"
    try:
        with RunLog("simZ", "boom", path=p):
            raise ValueError("bad")
    except ValueError:
        pass
    rows = list(csv.DictReader(p.open(encoding="utf-8")))
    assert rows[0]["outcome"].startswith("FAILED: ValueError")
