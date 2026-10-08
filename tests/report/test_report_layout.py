"""The PDF report keeps every table inside the card that holds it (seam: render_genome_report + WeasyPrint layout)."""
from __future__ import annotations

from pathlib import Path

import pytest
import yaml

weasyprint = pytest.importorskip("weasyprint")

from MATPredict.report.html import render_genome_report  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures"
CASES = sorted(p.stem for p in FIXTURES.glob("*.yaml"))
TOLERANCE_PX = 1.0   # sub-pixel rounding only


def _tables_outside_cards(html: str) -> list[tuple[str, float]]:
    """(table class, overshoot in CSS px) for every table whose border box ends right of its card's border box."""
    doc = weasyprint.HTML(string=html).render()
    found: list[tuple[str, float]] = []

    def walk(box, cards):
        el = getattr(box, "element", None)
        classes = (el.get("class") or "").split() if el is not None and hasattr(el, "get") else []
        if "card" in classes:
            cards = cards + [box]
        if getattr(box, "element_tag", None) == "table" and cards and not getattr(box, "is_table_wrapper", False):
            card = cards[-1]
            over = (box.border_box_x() + box.border_width()) - (card.border_box_x() + card.border_width())
            if over > TOLERANCE_PX:
                found.append((" ".join(classes), over))
        for child in getattr(box, "children", ()):
            walk(child, cards)

    for page in doc.pages:
        walk(page._page_box, [])
    return found


@pytest.mark.parametrize("case", CASES)
def test_no_table_runs_past_its_card_in_the_pdf(case):
    doc = yaml.safe_load((FIXTURES / f"{case}.yaml").read_text())
    html = render_genome_report(doc, sample=case, print_mode=True)
    assert _tables_outside_cards(html) == []
