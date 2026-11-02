"""Abstract and title handling for paper/manuscript.md and paper/PLAN.md.

The abstract has a 250-word limit, is reflowed to 100 columns and is copied into PLAN.md after
every revision; these steps were done by hand in four review rounds.
"""
import textwrap
from pathlib import Path

LIMIT = 250
WIDTH = 100
_M_HEAD, _M_TAIL = "## Abstract\n\n", "\n\n**Keywords"
_P_HEAD, _P_TAIL = "## Abstract (draft, 250 words max)\n\n", "\n\n## Section outline"


def _split(text, head, tail):
    before, rest = text.split(head, 1)
    body, after = rest.split(tail, 1)
    return before, body, after


def abstract(text):
    """The abstract of a manuscript as one line of words."""
    return " ".join(_split(text, _M_HEAD, _M_TAIL)[1].split())


def words(s):
    """Word count as journals count it: whitespace-separated tokens."""
    return len(s.split())


def reflow(s, width=WIDTH):
    return textwrap.fill(" ".join(s.split()), width, break_on_hyphens=False, break_long_words=False)


def set_abstract(text, new, limit=LIMIT):
    """Replace the manuscript abstract with `new`, reflowed; refuses to exceed the word limit."""
    n = words(new)
    if n > limit:
        raise ValueError(f"abstract has {n} words, limit {limit}")
    before, _, after = _split(text, _M_HEAD, _M_TAIL)
    return before + _M_HEAD + reflow(new) + _M_TAIL + after


def sync_plan(manuscript_md, plan_md):
    """Copy the manuscript's title and abstract into PLAN.md text; returns the new PLAN text."""
    title = manuscript_md.split("\n", 1)[0].lstrip("# ").strip()
    before, _, after = _split(plan_md, _P_HEAD, _P_TAIL)
    plan = before + _P_HEAD + reflow(abstract(manuscript_md)) + _P_TAIL + after
    i = plan.index("**Working title.**")
    j = plan.index("\n\n", i)
    return plan[:i] + "**Working title.** " + title + plan[j:]


def sync_files(root):
    root = Path(root)
    m = (root / "paper" / "manuscript.md").read_text(encoding="utf-8")
    p = root / "paper" / "PLAN.md"
    p.write_text(sync_plan(m, p.read_text(encoding="utf-8")), encoding="utf-8")
    return words(abstract(m))
