from pathlib import Path

import pytest

from thznerve.manuscript import LIMIT, abstract, set_abstract, sync_plan, words

M = "# A title\n\nAuthor\n\n## Abstract\n\nold words here\n\n**Keywords:** x\n\n## 1. Intro\n"
P = ("# Paper plan\n\n**Working title.** Old\ntitle\n\n**Target journal.** J\n\n"
     "## Abstract (draft, 250 words max)\n\nstale\n\n## Section outline\n\n1. x\n")


def test_set_abstract_reflows_and_keeps_the_rest():
    new = " ".join(["word"] * 60)
    out = set_abstract(M, new)
    assert abstract(out) == new
    assert all(len(line) <= 100 for line in out.split("\n"))
    assert out.startswith("# A title") and out.endswith("## 1. Intro\n")


def test_set_abstract_refuses_over_limit():
    with pytest.raises(ValueError):
        set_abstract(M, " ".join(["w"] * (LIMIT + 1)))


def test_sync_plan_copies_title_and_abstract():
    out = sync_plan(set_abstract(M, "new abstract"), P)
    assert "**Working title.** A title\n\n**Target journal.**" in out
    assert "## Abstract (draft, 250 words max)\n\nnew abstract\n\n## Section outline" in out


def test_real_manuscript_within_limit():
    m = (Path(__file__).resolve().parents[1] / "paper" / "manuscript.md").read_text(encoding="utf-8")
    assert words(abstract(m)) <= LIMIT
