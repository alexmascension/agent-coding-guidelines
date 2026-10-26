# Manuscript guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 1.0.0 (2026-10-26)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Each rule is followed by its rationale: ***Why***
> (mechanism), ***Trade-off*** (cost) and ***Deviate when*** (justified exceptions), where they apply.
> A rule without a key word is a default with the force of SHOULD. Rules carry stable identifiers (e.g.
> `MS-TRACE-001`); overrides refer to them (`CODING_GUIDELINES.md` §0, rule identifiers).

## 0. How to use this document

* Audience: agents and people who write, check and revise a manuscript (paper, preprint, thesis chapter, technical
  report) built on analyses run under these guidelines.
* Copy with `CODING_GUIDELINES.md` and `AGENTIC_GUIDELINES.md` into `guidelines/` of a project that produces a
  manuscript (`CODE-USE-001`). This guide applies the general rules on evidence, claims and verification of
  `AGENTIC_GUIDELINES.md` §3 to §6 to a manuscript and does not repeat them. What the Methods must report is in
  `STATISTICS_GUIDELINES.md` §12 (what to report) and in the reporting section of each domain guide used; figures
  will follow `FIGURE_GUIDELINES.md` (planned).
* Scope: the author's own manuscripts, and manuscripts whose authors have agreed to the procedure. A manuscript
  received for confidential peer review MUST NOT be processed by an AI service (`AGT-PRIV-003`), so §5
  (pre-submission review) does not apply to it.
* Sources: several rules adapt material of the Co-Scientist toolkit (Shocklab, University of Cape Town;
  <https://github.com/shocklab/co-scientist-toolkit>; text licensed under CC BY 4.0): `guidance/paper-CLAUDE.md`
  (itself adapted from Dominik Lukeš's agents-for-reproducibility course), the *Writing* and *Research papers*
  sections of `guidance/CLAUDE-starter.md`, and the skills `paper-review` (whose review templates are adapted from
  Dominik Lukeš's paper-reviewer skill, MIT), `claim-fact-checker` (adapted from Dominik Lukeš's agent skills,
  MIT), `reference-verifier` and `convergence-check`. They were rewritten as rules with identifiers, generalised
  beyond one agent and platform, and linked to the other guides of this repository.

## 1. Layout and inputs

* `MS-LAY-001` **The manuscript lives in the project repository**, in `paper/` (the default layout of
  `CODE-LANG-008` plus this folder), in a text format under version control: LaTeX, Markdown or Quarto.
  Bibliography files, the generated values file (`MS-TRACE-002`) and the claims table (`MS-CLAIM-001`) are in
  `paper/` too.
  * *Why:* a text source can be diffed, reviewed line by line and rebuilt by a script; changes to it are attributable
    and reversible like changes to code.
  * *Deviate when:* co-authors work in a word processor. Keep the exchanged files in `paper/` with a version in the
    name, record in the decisions log which file is current, and bring the numbers and figures in from the generated
    files (§2, traceability of results), never by retyping.
* `MS-LAY-002` **Papers and other documents the text relies on are read-only inputs**, kept in
  `data/literature/<citation key>/`: the original as received, when its licence allows a local copy, and the text
  converted from it on the machine (`AGT-EVID-006`, reading converted text). They are never edited
  (`CODE-PROV-001`, raw inputs) and never committed when their licence does not allow redistribution
  (`CODE-PROV-002`).
  * *Why:* a claim about a cited paper is checked against a fixed copy; a converted text kept with the original is
    read and searched in every later session without converting it again.
* `MS-LAY-003` **The terms of the manuscript are the controlled vocabulary of the project**
  (`PROJECT_CONTEXT.md` §4, domain vocabulary; `CODE-NAME-006`): a group, condition or cell type has the same name
  in the text, the figures, the tables and the code, or the vocabulary records the mapping.

## 2. Traceability of results

* `MS-TRACE-001` **Every figure and table is produced by a script from an explicit run**, never from "the latest"
  results and never edited by hand. The script takes the run directory as an argument and writes, next to the
  output, a sidecar `<name>.source.txt` naming the script, the run directory (`CODE-PROV-003`, run provenance) and
  the commit.
  * *Why:* two runs of the same script on different results produce figures that look the same; without the
    sidecar no one can say which run a submitted figure shows.
  * *Deviate when:* the final multi-panel layout is composed by hand. The scripted panels remain the source; the
    composite records the panels it contains and is composed again whenever a panel changes.
* `MS-TRACE-002` **Every number in the text comes from a generated file.** The analysis writes the values the text
  reports to a file the manuscript includes (LaTeX macros such as `\newcommand{\nDonors}{48}`, a Quarto or YAML
  variables file), rounded there and only there. A number that cannot be generated carries a comment in the source
  naming the file and run it was read from.
  * *Why:* numbers retyped into the text drift from the results with every rerun, and a rerun gives no error when
    the text keeps the old value (`AGT-CLAIM-003`, regenerated numbers).
* `MS-TRACE-003` **Broken or superseded runs are moved out of reach**: to `results/_ARCHIVE/<run>/`, with a
  `DO_NOT_USE.md` stating why and which run replaces them. Figure and table scripts MUST fail when given an archived
  run.
  * *Why:* a superseded run left in place is eventually read by a script whose argument was not updated, and the
    figure it produces is indistinguishable from a current one.
* `MS-TRACE-004` **Before each submission, regenerate every figure, table and number from a clean checkout** of the
  commit to be submitted, compare them with the manuscript, and tag the commit (e.g. `submission-<venue>-v1`); the
  tag and the runs it uses go in the decisions log.
  * *Why:* a manuscript assembled over months mixes outputs of different code versions; regenerating from one commit
    shows the mismatches before reviewers do.
* `MS-TRACE-005` **Data and code availability statements name what a reader needs to obtain them**: accessions for
  every dataset used (with the release), and the code repository with the submitted tag or an archived release with
  a persistent identifier (e.g. a DOI). The statements match the sources recorded under `CODE-PROV-002` (how to
  obtain inputs).

## 3. Claims

* `MS-CLAIM-001` **Keep a claims table for the load-bearing claims** of the title, abstract, contribution list,
  discussion and conclusions, in `paper/CLAIMS.md`: the claim as written, where in the manuscript it is made, the
  result that supports it (figure, table, section, generated value) and the kind of evidence it provides
  (`AGT-CLAIM-004`, kinds of evidence). A claim without a supporting result is rewritten, supported or removed.
  * *Why:* the claims a paper is cited for are in its abstract and conclusions, which are written last and drift
    furthest from the results; a table that ties each claim to its result makes the drift visible.
  * *Trade-off:* the table must be maintained through revisions; it covers load-bearing claims only, not every
    sentence.
* `MS-CLAIM-002` **The scope of each claim matches the design that produced it**: the population, tissue, platform
  and conditions studied. Causal language only where the design licenses it; generalisation beyond the studied
  population stated as a hypothesis; "first", "only", "novel", "establishes" and "X-fold" justified with what was
  checked (`AGT-CLAIM-001`, bounded claims).
* `MS-CLAIM-003` **Exploratory analyses are labelled as exploratory in the text**, and departures from the planned
  analysis (`STATISTICS_GUIDELINES.md` §1, design before analysis; `AGT-LOOP-005`, predictions fixed before the
  pilot) are reported as departures, with the reason. In the decisions log, exploratory method changes, departures
  from a fixed protocol and repaired tests (`AGT-VER-003`) are marked as such.
* `MS-CLAIM-004` **The limitations section states the known limitations**, including those recorded during the
  work (confounding that could not be avoided, `STAT-DES-004`; controls retired; populations not covered), and the
  conclusions do not contradict it.
* `MS-CLAIM-005` **The Methods report what `STATISTICS_GUIDELINES.md` §12 (what to report) and the reporting
  section of each domain guide require**, and the community reporting standard of the study type where one exists
  (e.g. ARRIVE for animal research [Percie du Sert 2020†], minSCe for single-cell RNA-seq [Füllgrabe 2020], PRISMA
  for systematic reviews and meta-analyses [Page 2021]). Look up the version the venue requires.

## 4. Citations

* `MS-CITE-001` **MUST NOT add a reference from memory or from model output without retrieving its record.** Every
  reference comes from a bibliographic record (DOI, PubMed, arXiv) or the publisher's page, and its existence and
  metadata (authors, title, year, venue, volume, pages, DOI) are checked against that record.
  * *Why:* language models produce references that do not exist and introduce errors into references that do: in
    short literature reviews generated by ChatGPT, 18% of GPT-4's references were fabricated and 24% of the real
    ones had substantive errors [Walters 2023†].
* `MS-CITE-002` **Before submission, check the whole bibliography against bibliographic databases** (e.g. Crossref,
  OpenAlex, Semantic Scholar, arXiv), with a tool rather than by reading, and report: references not found,
  metadata that disagree with the record, DOIs that do not resolve, works cited under two keys, and, from the
  source, keys cited but not listed and listed but not cited. When a reference list must first be transcribed
  from a PDF or a word processor file, copy each entry as printed, without correcting it, since the check is meant
  to find its errors.
* `MS-CITE-003` **Each citation supports the sentence it is attached to.** A reference that exists is not thereby
  support. For load-bearing claims, read the cited passage in the cited work's text (not its abstract;
  `AGT-EVID-004`, primary sources) and classify the citation: supported, supported only through a secondary source,
  unsupported, contradicted, misattributed (the cited work does not say it), or a quotation that is not verbatim.
  Check long lists one claim per fresh context (`AGT-DEL-005`).
  * *Why:* in a meta-analysis of medical journal articles, about a quarter of quotations of cited work contained an
    error, and about one in eight a major one [Jergas 2015†].
  * *Trade-off:* reading the cited passage for every claim is expensive; scale it to the load-bearing claims of
    `MS-CLAIM-001` (claims table).
* `MS-CITE-004` **Cite the primary source of a finding**, not a review or a later paper that reports it; when only
  a secondary source can be cited, say so in the text. Cite datasets by accession and by the paper that describes
  them, and software with its version and its requested citation.
* `MS-CITE-005` **Read cited papers within the limits of `AGENTIC_GUIDELINES.md` §10** (confidential and licensed
  material): open copies in full; for subscription papers without an open copy, the title, abstract and short
  quotes, unless the institution's agreement permits more (`AGT-PRIV-004`).

## 5. Pre-submission review

* `MS-REV-001` **Before submission, review the manuscript along independent dimensions**, each by a reviewer with a
  fresh context (`AGT-CHAL-002`, fresh-context review):
  * Methods and statistics: whether the design tests the claims, controls and baselines (fairly tuned, not
    strawmen), experimental units and replicates, multiple testing, leakage, reproducibility from the paper alone,
    and figures whose scales or selection misrepresent the data.
  * Writing: structure, whether the contribution is clear early, terms explained before use, redundancy.
  * Positioning: whether cited work is represented accurately, missing related or concurrent work, and novelty
    claims checked against the retrieved text of prior work, not recalled.
  * Internal consistency: abstract against results, contributions promised against delivered, numbers in the text
    against tables and figures, captions against content, notation, cross-references.
  * The strongest objection to the main claim, with what evidence or argument would defeat it.
  * Proofs and derivations, line by line, where the manuscript has them.
* `MS-REV-002` **Each finding states its location** (section, paragraph, figure, table, equation), severity
  (critical, major, minor, style), confidence and a suggested remedy. Findings are deduplicated and ranked by
  severity, then confidence, then dependency (a fix that makes another unnecessary comes first). Disagreements
  between reviewers go to the author with both sides, not averaged.
* `MS-REV-003` **Review against the venue's current criteria**: look up its reviewer guidelines and author
  instructions (dimensions and scales, word limits, required checklists, data and code policy, policy on AI use)
  rather than recalling them (`AGT-EVID-003`, look-ups), and state when generic criteria had to be used instead.
* `MS-REV-004` **MUST NOT apply review changes without the author's approval of each item.** Approved changes are
  applied as the smallest edit that addresses the item, in the author's voice; nothing outside the approved list
  is changed.
  * *Why:* the author is responsible for the text and is the domain expert; a reviewer's suggestion applied
    unapproved can be wrong in ways only the author sees.
* `MS-REV-005` **After applying changes, verify them**: every approved change is implemented as specified, no
  unapproved change was introduced, no new inconsistency appeared (a claim changed in one section and not in
  another), and the manuscript builds and is checked as rendered, not only in source (`AGT-VER-001`, outcome over
  proxy).

## 6. Revisions and responses

* `MS-RESP-001` **A question from a reviewer, editor or co-author is answered, not acted on** (`AGT-USER-001`): the
  author decides which points lead to changes.
* `MS-RESP-002` **The response letter quotes each point, answers it and gives the location of each change** (section,
  page and line, or the changed text). A point that is not addressed says why.
* `MS-RESP-003` **Each revision is a tagged commit** (`MS-TRACE-004`, submission tags), so that a marked-up
  difference between submitted versions can be produced (e.g. with `latexdiff`) and every figure and number of each
  version can be regenerated.
* `MS-RESP-004` **Check whether revision rounds converge** (`AGT-LOOP-002`, convergence of revision loops): when
  each round raises a new kind of problem, or the abstract has been reframed three or more times on unchanged
  results, say that the bottleneck is the data or the question, not the wording, and propose a narrower claim, more
  evidence or a different question.

## 7. Writing

* `MS-WRITE-001` **Say what is computed and what it gives.** Avoid figurative verbs ("pushes the error down", "sits
  inside the bound"); explain a term the reader may not know on first use and define each acronym once; give
  numbers with their units and their universe (`CODE-NUM-009`).
* `MS-WRITE-002` **Write for the reader who acts on the text.** Documents for funders, committees and new
  collaborators state the facts and the current plan (what happens, when, and who does it), without persuasion,
  reassurance or accounts of how the plan changed; open questions are kept.
* `MS-WRITE-003` **Text drafted or revised by an AI is the authors' responsibility**, checked by them like any
  other draft, and its use is disclosed as the venue requires (`MS-REV-003`, venue criteria). AI tools are not
  authors.
* `MS-WRITE-004` **A rewrite by a second model or service gets a brief and is checked against it.** The brief
  states the reader, the aim and the choices that must survive; the rewrite is compared with both the draft and the
  brief before the author sees it. Only material that `AGENTIC_GUIDELINES.md` §10 (confidential and licensed
  material) allows to leave is sent.
* `MS-WRITE-005` **Give an honest assessment, not compliance.** When asked to improve a manuscript, say whether the
  work answers an interesting question, not only whether it is publishable, and push back on a change that weakens
  the evidence or overstates it.

## 8. What to record in `PROJECT_CONTEXT.md`

1. `MS-REP-001` The target venue, the date its criteria and author instructions were consulted, and its policy on
   AI use (§5, pre-submission review).
2. `MS-REP-002` The manuscript format and build command, the generated values file, the claims table and the
   location of figure and table scripts (§1, layout and inputs; §2, traceability of results).
3. `MS-REP-003` Submission and revision tags with the runs each uses (`MS-TRACE-004`, `MS-RESP-003`).
4. `MS-REP-004` How AI tools were used in writing and review, for the disclosure the venue requires (`MS-WRITE-003`).
