"""Supplement scan (layer 3a): unit splitting, gene-context filter, notation regexes."""

from __future__ import annotations

import io
import zipfile

import openpyxl

from tessera.supplements import scan, units


def _zip(files: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for k, v in files.items():
            z.writestr(k, v)
    return buf.getvalue()


def _xlsx(rows) -> bytes:
    wb = openpyxl.Workbook()
    for r in rows:
        wb.active.append(r)
    b = io.BytesIO()
    wb.save(b)
    return b.getvalue()


def test_row_level_gene_context_including_spaced_symbol():
    # row text as it appears in a real supplement: "SETD 1A | c.4596_4597insG | c.4582-2delAG"
    z = _zip({"s2.xlsx": _xlsx([["Gene", "Mutation 1", "Mutation 2"],
                               ["SETD 1A", "c.4596_4597insG ", "c.4582-2delAG "],
                               ["OTHERGENE", "c.100A>G", "p.Arg34Ter"]])})
    hits = scan("1", units(z), ["SETD1A"])
    assert {h.notation for h in hits} == {"c.4596_4597insG", "c.4582-2delAG"}
    assert all(h.where == "sheet=Sheet;row=2" for h in hits)


def test_protein_frameshift_shorthand_is_not_a_hit_but_substitution_is():
    z = _zip({"t.csv": b"gene,change\nSETD1A,p.Leu1533fs\nSETD1A,p.Arg990Ter\nSETD1A,p.(Gln737*)\n"})
    notations = {h.notation for h in scan("1", units(z), ["SETD1A"])}
    assert notations == {"p.Arg990Ter", "p.Gln737*"}


def test_transcript_named_in_row_is_kept():
    z = _zip({"t.csv": b"SETD1A,NM_014712.3,c.2968C>T\n"})
    (h,) = scan("1", units(z), ["SETD1A"])
    assert h.transcript == "NM_014712.3" and h.notation == "c.2968C>T"


def test_unreadable_member_is_recorded_not_fatal():
    us = units(_zip({"broken.xlsx": b"not really xlsx"}))
    assert us and us[0].where == "unreadable"


def test_other_genes_variants_on_the_same_page_are_not_attributed():
    # shape of a real PDF supplement page: several genes' variants in one text block
    page = ("Refractory SETD1A c.4582-2delAG p.? de novo LP 1056 Complex febrile N/A SCN8A c.3955G>T "
            "p.Ala1319Ser CADD (23.7) SCN3A c.5147C>T p.Val1532Leu REVEL (0.7) PVS1 PM2")
    z = _zip({"t.csv": ("x," + page + "\n").encode()})
    assert {h.notation for h in scan("1", units(z), ["SETD1A"])} == {"c.4582-2delAG"}
