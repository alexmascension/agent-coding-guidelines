"""Tests of the BibTeX checks of ``tools/check_guides.py``."""

from tools.check_guides import bib_problems

BIB = (
    "@article{Smith_2020, title={A}, author={Smith, Ann and Lee, Bo}, year={2020}, month=sep, "
    "DOI={10.1/ABC} }\n"
    "@article{Jones_1986, title={B}, author={{Some Consortium} and Jones, Cy}, year={1986}, "
    "month={Oct}, DOI={10.1016/0197-2456(86)90046-2} }\n"
)
MD = (
    "- **[Smith 2020]** Smith & Lee (2020). A. *J*. [doi:10.1/abc](https://doi.org/10.1/abc)  \n"
    "  Note.\n"
    "- **[Jones 1986]** Jones (1986). B. *J*. "
    "[doi:10.1016/0197-2456(86)90046-2](https://doi.org/10.1016/0197-2456(86)90046-2)  \n"
)


def test_consistent_files_have_no_problems():
    assert bib_problems("X", BIB, MD) == []


def test_reports_key_used_twice():
    bib = BIB.replace("Jones_1986", "Smith_2020")
    assert bib_problems("X", bib, MD) == ["X.bib: key Smith_2020 used more than once"]


def test_reports_doi_listed_in_markdown_only():
    md = MD + "- **[New 2021]** New (2021). C. *J*. [doi:10.2/new](https://doi.org/10.2/new)  \n"
    assert bib_problems("X", BIB, md) == ["X.md: DOI 10.2/new has no entry in X.bib"]


def test_reports_doi_listed_in_bibtex_only():
    md = MD.split("- **[Jones")[0]
    assert bib_problems("X", BIB, md) == [
        "X.bib: DOI 10.1016/0197-2456(86)90046-2 has no entry in X.md"
    ]


def test_reports_doi_in_two_bibtex_entries():
    bib = (
        BIB
        + "@article{Smith_2020b, title={A}, author={Smith, Ann}, year={2020}, DOI={10.1/abc} }\n"
    )
    assert bib_problems("X", bib, MD) == ["X.bib: DOI 10.1/abc in more than one entry"]


def test_reports_empty_author_left_by_group_author():
    bib = BIB.replace("author={Smith", "author={ and Smith")
    assert bib_problems("X", bib, MD) == ["X.bib:1: empty author name"]


def test_reports_month_that_is_not_a_macro():
    bib = BIB.replace("month=sep", "month=Sept")
    assert bib_problems("X", bib, MD) == ["X.bib:1: month=Sept is not a BibTeX month macro"]
