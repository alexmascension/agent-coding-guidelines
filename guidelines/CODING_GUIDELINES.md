# Coding guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 3.4.0 (2026-10-24)
>
> The key words **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT** and **MAY** are used as in
> RFC 2119. MUST rules are hard constraints, reserved for conditions of validity that can be checked. SHOULD rules
> are strong defaults (methods, thresholds, tools); deviating from them requires a reason stated in the change or in
> `PROJECT_CONTEXT.md`. When a MUST rule is followed by *Deviate when*, the situation described there is outside the
> scope of the rule, not an exception to it. A rule stated without a key word (e.g. "**Default: Harmony**") is a
> default with the force of SHOULD.
>
> Every rule has a stable identifier (e.g. `SC-QC-004`) written after its list marker; overrides and
> discussions refer to rules by identifier (§0, rule identifiers).

## 0. How to use this document

* `CODE-USE-001` This file is project-agnostic. Copy it **verbatim** into the folder `guidelines/` at the root of every
  project, together with the domain guidelines and, for clusters, the site profile (`guidelines/hpc_sites/`); the
  reference files the guides cite go to `guidelines/references/`. In the upstream repository, guides and site profiles
  are under `guidelines/` and reference files under `references/`. `AGENTS.md`, `PROJECT_CONTEXT.md` and `README.md`
  stay in the project root.
* `CODE-USE-002` Everything specific to a project (stack details, layout, vocabulary, data invariants,
  conventions, interfaces, decisions) goes into `PROJECT_CONTEXT.md` (see §2, documentation
  set).
* `CODE-USE-003` Do not edit the copy of this file to adapt it to a project. Record adaptations in the
  *Overrides* section of `PROJECT_CONTEXT.md`. This keeps the copy comparable with upstream and
  lets it be updated without losing project decisions.
* `CODE-USE-004` When this file and `PROJECT_CONTEXT.md` conflict, `PROJECT_CONTEXT.md` wins **only** for
  rules it explicitly lists under *Overrides*. Otherwise this file wins.
* `CODE-USE-005` Domain guidelines from the same repository (statistics, figures, RNA-seq, …) are copied
  alongside this file when relevant to the project and follow the same override mechanism.
* `CODE-USE-006` **Precedence between guidelines**: `PROJECT_CONTEXT.md` *Overrides*, then the most specific domain guideline for
  the question at hand (e.g. `SEX_INFERENCE_GUIDELINES.md` on reference routing over `RNASEQ_BULK_GUIDELINES.md`),
  then the more general one, then this file. A conflict between domain guidelines that this order does not resolve
  is decided in *Overrides* and reported upstream. The adopted guidelines and their versions are listed at the start
  of *Overrides* (§2, required structure).
* `CODE-USE-007` **Versioning.** Every guideline file carries a semantic version, and changes are recorded in
  the repository `CHANGELOG.md`:
  * MAJOR: a rule is removed (its identifier retired) or reversed (what it requires is inverted or substantially
    changed), so existing overrides may no longer apply.
  * MINOR: new rules, or sections reorganized or renumbered without changing rules (identifiers are kept).
  * PATCH: wording or clarifications that do not change what a rule requires.
* `CODE-USE-008` **Rule identifiers.** Every rule (top-level list item of a section with rules) carries a stable
  identifier `PREFIX-TOPIC-NNN` written as inline code after its list marker, e.g. `SC-QC-004`: the prefix names
  the guide (`CODE`, `STAT`, `BULK`, `SC`, `SP`, `META`, `SEX`, `ID`, `DECON`, `HPC`), the topic names the section where the
  rule was created, and the number is sequential within the topic.
  * Identifiers are never renumbered or reused. A rule moved to another section keeps its identifier; a removed
    rule's identifier is marked `retired` in `RULE_IDS.tsv`; a new rule takes the next free number of its topic
    (`tools/rule_ids.py assign`).
  * Overrides, decisions and discussions refer to rules by identifier (the section number MAY be added for
    readability: `SC-QC-003` (§4.2)).
  * *Why:* section numbers change when text is added or reorganized; an override that points to "§4.2, second
    bullet" silently points to another rule after an edit, while an identifier does not move.

## 1. Priority order

When rules or goals conflict, prioritize:

**correctness > readability > backward compatibility > consistency with the existing codebase >
brevity > performance**

## 2. Project documentation set

Every project developed with agents MUST contain at least:

| File | Content |
|---|---|
| `AGENTS.md` | Entry point read automatically by coding agents. Pointers only (template below). |
| `PROJECT_CONTEXT.md` | Project-specific decisions and constraints for development (structure below). |
| `guidelines/CODING_GUIDELINES.md` | This file, unmodified (with the copied domain guidelines and references, §0). |
| `README.md` | What the project does, how to install it, how to run it, where outputs go. |

Add further documents only when they carry information that does not fit the above, e.g.
`PROJECT_SCOPE.md` (scientific scope and design rationale), `DECISIONS.md` (if the decisions log
outgrows `PROJECT_CONTEXT.md`), copied domain guidelines in `guidelines/` (`STATISTICS_GUIDELINES.md`,
`FIGURE_GUIDELINES.md`, `RNASEQ_BULK_GUIDELINES.md`, …), and, for repositories that are published,
`LICENSE`, `CITATION.cff` and `CHANGELOG.md`.

### `AGENTS.md` template

`AGENTS.md` MUST NOT duplicate rules or project facts; it only points to the files that hold
them, so there is a single source of truth.

```markdown
# Agent instructions

Before any work in this repository, read and follow, in this order:

1. `PROJECT_CONTEXT.md` — project-specific decisions; its *Overrides* section takes
   precedence over the guidelines for the rules it lists.
2. `guidelines/CODING_GUIDELINES.md` — general development rules.
3. Domain guidelines in `guidelines/` (`*_GUIDELINES.md`; site profiles in `guidelines/hpc_sites/`), when the task
   touches them.
```

For agents that do not read `AGENTS.md` (e.g. Claude Code reads `CLAUDE.md`), add the file the
agent expects containing only a reference to `AGENTS.md` (in `CLAUDE.md`: `@AGENTS.md`).

### Required structure of `PROJECT_CONTEXT.md`

`PROJECT_CONTEXT.md` MUST keep the following sections, in this order. A section that does not
apply MUST be kept with the text `None.` so its absence is an explicit decision rather than an
omission.

1. `CODE-PROJ-001` **Purpose** — one paragraph; link to `README.md` / `PROJECT_SCOPE.md` for detail.
2. `CODE-PROJ-002` **Stack** — languages, key libraries, external tools **with versions**, environment name and
   how to create it.
3. `CODE-PROJ-003` **Layout** — directories, what goes in each, where inputs and outputs live, which paths
   these guidelines apply to and which are exempt.
4. `CODE-PROJ-004` **Domain vocabulary** — controlled labels, categorical values and column names, and the
   module where they are declared (see §6, controlled vocabulary).
5. `CODE-PROJ-005` **Data model and invariants** — atomic tables, entity keys, derivation rules, filters that
   define each table (see §7, tables and schemas).
6. `CODE-PROJ-006` **Conventions** — coordinate systems, genome build and annotation release, units,
   normalization, thresholds and their rationale (see §11, numerical code).
7. `CODE-PROJ-007` **Interfaces** — outputs consumed by notebooks, other repositories or people; sibling
   libraries that MUST be reused rather than copied (see §21, interfaces).
8. `CODE-PROJ-008` **Overrides** — first, the version of each adopted guideline file; then every deviation
   from them, each with the rule identifier (e.g. `SC-QC-003`), the replacement rule and the reason.
9. `CODE-PROJ-009` **Decisions log** — dated entries: decision, rationale, affected outputs. MAY instead link to
   `DECISIONS.md`.

Keep `PROJECT_CONTEXT.md` synchronized with the code. A change that alters a convention, an
interface, a threshold or a data invariant MUST update it in the same change.

## 3. Languages and tooling

* `CODE-LANG-001` Python 3.11+ is the default language.
* `CODE-LANG-002` Prefer a maintained Python implementation of an established method over calling R or another
  language. Use it only if it is validated against the reference implementation; when results
  must match the reference (published method, comparison with prior work), state which
  implementation and version was used, and document known differences.
  * *Why:* one language keeps the environment, the data structures and the tests in one place, and Python code is
    easier for the project's readers to inspect. *Trade-off:* ports can differ numerically from the reference;
    hence the validation requirement. *Deviate when:* no port exists, or the port is not validated for the analysis
    at hand.
* `CODE-LANG-003` Use R, bash or other languages only when no adequate Python option exists or when the user
  requests it. Isolate such code in standalone scripts with a file-based interface (inputs and
  outputs on disk, parameters on the command line) called from Python. Avoid in-process bridges
  such as `rpy2` unless there is a concrete reason.
  * *Why:* a file interface makes every intermediate result inspectable and lets each language run in its own
    environment; in-process bridges tie the versions of two runtimes together and fail in ways that are hard to
    reproduce.
* `CODE-LANG-004` Bash is for glue only (calling tools, moving files, little data manipulation):
  * Rewrite a script in Python once it exceeds about 100 lines or needs non-trivial control
    flow.
  * Start every script with `set -euo pipefail`; quote all variable expansions.
  * Use `$(…)` rather than backticks and `[[ … ]]` rather than `[ … ]`; send errors to STDERR;
    put the logic in a `main` function called at the end of the file.
  * Scripts SHOULD pass `shellcheck`.
  * Add `set -E` and an `ERR` trap that reports the line, the command and the exit code to STDERR:
    `trap 'echo "ERROR line ${LINENO}: ${BASH_COMMAND} (exit $?)" >&2' ERR`. Do not discard STDERR (`2>/dev/null`)
    of a command whose failure aborts the script; tolerate expected failures explicitly (`|| true`).
    * *Why:* under `set -e`, a failing command whose STDERR is discarded ends the script with a non-zero code and an
      empty log. The trap does not fire in `if`/`while` conditions or `&&`/`||` lists; for a failing pipeline it
      names the last command of the pipeline, not necessarily the one that failed.
* `CODE-LANG-005` R code, when needed, follows the tidyverse style guide, is formatted with `styler` and linted
  with `lintr`.
* `CODE-LANG-006` Default Python stack (a project MAY override it in `PROJECT_CONTEXT.md`): polars for tables,
  NumPy/SciPy for numerics, `statsmodels` for statistical models, matplotlib/seaborn for figures,
  pytest for tests, ruff for linting and formatting.
* `CODE-LANG-007` Environment: a conda environment plus `pip install -e .` of the project package. Export the
  resolved environment (`conda env export --no-builds > environment.lock.yml` or `conda-lock`)
  whenever dependencies change.
  * *Why:* an environment file with version ranges resolves to different packages months later; the lockfile
    records what actually ran.
* `CODE-LANG-008` Default layout (a project MAY override it):
  `src/<package>/` library code, `src/scripts/` command-line drivers, `src/notebooks/`
  analysis notebooks, `tests/`, `data/` inputs, `results/` outputs.
  These guidelines apply to library code, scripts and tests; notebooks follow §22 (notebooks).

## 4. General principles

* `CODE-PRIN-001` Write clean, readable and maintainable code.
* `CODE-PRIN-002` Prefer simple, explicit solutions over clever or unnecessarily complex ones.
* `CODE-PRIN-003` Give each function one clear responsibility; split large or deeply nested logic into focused
  helpers when this improves clarity. A function longer than about 40 lines is a signal to
  consider splitting it, not a hard limit.
* `CODE-PRIN-004` Extract repeated logic into reusable helpers; do not introduce abstractions unless they
  genuinely simplify the code.
* `CODE-PRIN-005` Use descriptive names and intermediate variables when they make logic easier to follow.
  Avoid compact one-liners that reduce readability; a comprehension SHOULD have at most one
  `for` clause and one condition.

## 5. Changes to existing code

* `CODE-CHG-001` Preserve existing behaviour and interfaces unless the requested change requires otherwise.
* `CODE-CHG-002` Make the smallest coherent change that fully solves the requested problem.
  * *Why:* small diffs can be reviewed; unrelated changes mixed into a fix hide regressions and make reverting one
    change impossible without reverting the other.
* `CODE-CHG-003` MUST NOT refactor or rewrite unrelated code while implementing a feature or bug fix.
* `CODE-CHG-004` Reuse existing project utilities and conventions instead of creating parallel implementations.
* `CODE-CHG-005` If a broader refactor would materially improve the solution, propose it separately.
* `CODE-CHG-006` If you notice an unrelated bug, report it; do not fix it in the same change.

## 6. Naming and controlled vocabulary

* `CODE-NAME-001` `UPPER_SNAKE_CASE` for module-level and configuration constants; `lowercase_with_underscores`
  for functions, parameters, variables and other non-constant names.
* `CODE-NAME-002` Prefer descriptive, self-contained names over abbreviations, unless the abbreviation is
  standard in the domain (`tpm`, `fdr`, `cds`).
* `CODE-NAME-003` Loop and comprehension variables stay descriptive when their meaning is not obvious; avoid
  `x`, `a`, `tmp`.
* `CODE-NAME-004` Short names are allowed only for standard mathematical notation (`n`, `k`, `p`, `X`), for
  notation that follows a reference paper or algorithm (cite the source in the docstring), or
  for labels defined as such in `PROJECT_CONTEXT.md`. Public APIs keep descriptive names.
* `CODE-NAME-005` For hierarchical information, order name components from broader to more specific
  (`gene_copy_count`, not `count_copy_gene`).
* `CODE-NAME-006` **Controlled vocabulary has a single source of truth.** Categorical labels, region names,
  classification values and column names MUST be declared once (a `config.py` or `schema.py`
  module, as listed in `PROJECT_CONTEXT.md` §4, domain vocabulary) and imported everywhere
  else. MUST NOT invent synonyms, aliases or alternative spellings in modules or notebooks.
  * *Why:* two spellings of one label (`XAR` and `X_added`) silently split one group into two, and every downstream
    count is wrong without any error being raised.

## 7. Data, tables and schemas

* `CODE-DATA-001` polars is the default for every table. pandas appears only at boundaries that require it
  (plotting libraries, `statsmodels`, third-party APIs): convert once, as late as possible, and
  never round-trip.
* `CODE-DATA-002` MUST NOT migrate existing code between pandas and polars solely for stylistic reasons.
* `CODE-DATA-003` Schemas are part of the function contract: document required columns and important dtypes
  (see §9, documentation), and validate them at public entry points.
* `CODE-DATA-004` Preserve row order, schemas and dtypes unless the operation explicitly changes them. In
  polars, `group_by` does not preserve group order; use `maintain_order=True` or sort explicitly
  when order matters downstream.
* `CODE-DATA-005` Joins MUST state their expected cardinality (polars `validate="1:1"` / `"m:1"`) when a
  violation would silently duplicate or drop rows, and SHOULD check row counts before and after.
  * *Why:* a many-to-many join on a key that was assumed unique duplicates rows without error, inflating every
    count computed afterwards.
* `CODE-DATA-006` Prefer columnar or vectorized operations for bulk transformations; do not force
  vectorization when row-oriented logic is clearer or required by an external API.
* `CODE-DATA-007` **Aggregates are derived, never accumulated.** Keep the most granular (atomic) table as the
  source of truth. Re-derive summaries from it after any filter. Combine previously aggregated rows only when the
  groups are disjoint and the statistic is decomposable (counts, sums; means only together with their
  denominators), and state the cardinality of each group.
  * *Why:* summing pre-aggregated rows double counts entities that appear in several groups and keeps results
    computed under old filters; re-deriving from the atomic table guarantees that every summary reflects the
    current filters.
* `CODE-DATA-008` **Entity keys are stable.** Each entity has a documented key (`PROJECT_CONTEXT.md` §5, data
  model) that is identical across modules. Any rule that merges or splits entities is a named
  constant and is applied in one place.
* `CODE-DATA-009` Use Parquet for intermediate and output tables (it preserves dtypes). Use TSV/CSV only for
  human-facing exports or tool interoperability.
  * *Why:* CSV loses dtypes: identifiers with leading zeros become numbers, nulls become empty strings, and
    categorical columns are re-parsed on every load.

## 8. Type hints

* `CODE-TYPE-001` Add type hints to public functions and non-trivial internal functions.
* `CODE-TYPE-002` Use built-in generics (`list[str]`, `dict[str, int]`) and `X | None`, not `typing.List` or
  `Optional[X]`.
* `CODE-TYPE-003` Do not add annotations that reduce readability or that cannot express the type meaningfully.

## 9. Documentation

### Language and docstrings

* `CODE-DOC-001` Write repository documentation, notebook prose, code comments, docstrings, log and error
  messages, and figure/table labels in **English**, regardless of the conversation language.
* `CODE-DOC-002` Write NumPy-style docstrings for public functions and non-trivial internal functions:

  ```text
  Parameters
  ----------
  counts : pl.DataFrame
      One row per gene and sample. Required columns: ``gene_id`` (Utf8),
      ``sample_id`` (Utf8), ``count`` (UInt32).

  Returns
  -------
  pl.DataFrame
      Description, including columns added or removed.

  Raises
  ------
  ValueError
      If a required column is missing.
  ```

* `CODE-DOC-003` For functions that take or return tables, the docstring MUST list the required columns and
  their dtypes, and state any ordering or uniqueness assumptions.
* `CODE-DOC-004` Document assumptions, invariants, edge cases and non-obvious behaviour.

### Self-contained comments and documents

* `CODE-DOC-005` Every comment and document MUST be understandable by a reader who knows the project scope but
  has no knowledge of how the code was developed.
* `CODE-DOC-006` Comments explain *why*, not *what*. Use them for non-obvious decisions, assumptions,
  performance considerations and tricky logic; do not comment self-explanatory code.
* `CODE-DOC-007` State the **reason**, not the **history**. MUST NOT refer to the conversation, request,
  meeting, previous version or decision that led to the code (`as discussed`, `per user
  request`, `changed after review`, `previously we used X`). Write the reason itself
  (`# Merge hits closer than MERGE_HITS_BP: aligners split one copy at short indels`).
  Decision history belongs only in the decisions log.
* `CODE-DOC-008` MUST NOT leave commented-out code; version control keeps old code.
* `CODE-DOC-009` `TODO` comments use the format `# TODO: <issue or link> - <description>` and are reported to
  the user when created.
* `CODE-DOC-010` Documents that define rules (guidelines, `PROJECT_CONTEXT.md`) give each non-obvious rule its rationale: the
  mechanism (*Why*), its cost (*Trade-off*) and the justified exceptions (*Deviate when*).
  * *Why:* a rule whose reason is visible can be checked and applied to cases it does not mention; a bare rule
    can only be obeyed or broken.
* `CODE-DOC-011` Keep documentation minimal and current: delete or update documents that no longer describe
  the code. Outdated documentation is worse than missing documentation.
* `CODE-DOC-012` Keep documentation synchronized with the implementation, including `README.md` and
  `PROJECT_CONTEXT.md` when an interface or design decision changes.

### Markdown

* `CODE-DOC-013` One H1 per file; ATX headings (`#`), each heading unique within the file.
* `CODE-DOC-014` Fenced code blocks always declare their language (`python`, `bash`, `text`).
* `CODE-DOC-015` Inline equations use `$...$`; display equations use `$$...$$` on their own lines.
  MUST NOT use `\(...\)` or `\[...\]`. This also applies to notebook Markdown cells.
* `CODE-DOC-016` **Cross-references state what they point to.** A reference to a section, rule, figure, table
  or numbered item MUST be followed by a short description, e.g. `§13 (data provenance)`,
  `Fig. 2 (copy-age distribution)`, `step 4 (normalization)`. The description need not be the
  exact title, only enough to identify the target without opening it.

## 10. Robustness

* `CODE-ROB-001` Validate important inputs and fail early with error messages that precisely match the failing
  condition, show interpolated values unambiguously (`f"Not a probability: {p=}"`) and are easy
  to grep.
* `CODE-ROB-002` MUST NOT use `assert` to validate inputs or enforce program logic: assertions are removed
  under `python -O`. Use `assert` only in tests and for internal invariants whose removal would
  not change behaviour.
* `CODE-ROB-003` Handle expected edge cases explicitly (empty inputs, missing groups, single-sample designs).
* `CODE-ROB-004` MUST NOT silently swallow errors. Avoid broad `except Exception` unless there is a specific,
  documented reason, and re-raise or log with context.
* `CODE-ROB-005` MUST NOT use mutable objects (`[]`, `{}`, a DataFrame) as default argument values; use `None`
  and create the object inside the function.
  * *Why:* Python evaluates default values once, at function definition, so a mutable default is shared by all
    calls and accumulates state between them.
* `CODE-ROB-006` MUST NOT mutate input arguments (DataFrames, arrays, lists, dicts) in place unless this is
  intentional and documented.
* `CODE-ROB-007` Open files and other resources with `with` statements.
* `CODE-ROB-008` Avoid "power features" (metaclasses, monkeypatching, dynamic attribute access via
  `getattr`/`setattr` on computed names, import-time side effects) unless there is no simpler
  alternative; document why when used.

## 11. Numerical and scientific code

* `CODE-NUM-001` Make units, coordinate systems, normalization conventions and statistical assumptions explicit
  when they are not obvious.
* `CODE-NUM-002` Avoid exact floating-point equality unless mathematically guaranteed; use tolerances.
* `CODE-NUM-003` Handle NaN, null, infinity and empty inputs deliberately. In polars, `null` and `NaN` are
  distinct; decide which one each column uses.
* `CODE-NUM-004` Check for `None` with `is None` / `is not None`, never by truthiness: `0`, `0.0` and empty
  containers are falsy, and the truth value of an array or DataFrame is ambiguous or raises.
* `CODE-NUM-005` MUST NOT silently change numerical precision or relevant dtypes.
* `CODE-NUM-006` Prefer established numerical, statistical and domain implementations over custom ones.
* `CODE-NUM-007` **One internal coordinate convention.** Use a single internal convention (default: 0-based,
  half-open). Convert every external format at the parser, at load time, never later. Record the
  convention of each input format in `PROJECT_CONTEXT.md` §6 (conventions).
  * *Why:* off-by-one errors arise at every conversion between 0- and 1-based systems; converting once, at the
    parser, confines the risk to one tested function.
* `CODE-NUM-008` Genome build and annotation release are part of every coordinate-bearing table and are
  checked when tables are combined.
* `CODE-NUM-009` **Every reported number states its universe**: which input set, which method or branch, and
  which filters it was computed on. A count without its universe is not a result.
* `CODE-NUM-010` **Classifications carry their rule.** A column produced by a documented classification rule
  (as opposed to a direct observation) MUST be reported together with the identifier or version
  of the rule that produced it.
  * *Why:* the same column name produced by two versions of a rule looks comparable but is not; without the rule
    version, results from different runs cannot be reconciled.
* `CODE-NUM-011` Thresholds and analysis parameters are named constants in one configuration module, each with
  a short comment giving its rationale or source. MUST NOT hard-code them inline.
  * *Why:* inline values get copied and then diverge between modules; a named constant with its rationale is
    changed in one place and is visible in review.
* `CODE-NUM-012` Statistical practice (test choice, multiple testing, effect sizes, reporting) follows
  `STATISTICS_GUIDELINES.md` when present.

## 12. External tools and subprocesses

* `CODE-EXT-001` Call external tools with `subprocess.run([...], check=True)` using an argument list.
  MUST NOT use `shell=True` with interpolated values.
  * *Why:* with `shell=True`, file names with spaces or special characters are split or interpreted by the shell,
    and interpolated values can execute arbitrary commands.
* `CODE-EXT-002` Check that required executables exist (`shutil.which`) before starting long work, and fail
  with a message naming the missing tool.
* `CODE-EXT-003` Log the exact command line and the tool version for every external call that produces a
  result; store the version with the outputs (see §13, provenance).
* `CODE-EXT-004` Capture stderr and include its tail in the error message when a tool fails.
* `CODE-EXT-005` Parse tool outputs in dedicated parser functions that convert to internal conventions (see
  §11, coordinate convention) and return typed tables.

## 13. Data management and provenance

* `CODE-PROV-001` Raw inputs (`data/raw/` or equivalent) are read-only. MUST NOT modify, overwrite or delete
  them.
* `CODE-PROV-002` MUST NOT commit data, results or large binaries to git. Record how to obtain each input
  (URL, accession, release, checksum) in `README.md` or a download script.
  * *Deviate when:* small derived tables (about 5 MB per file or less) needed to reproduce key results without
    re-running expensive upstream steps (validated run tables, sample manifests, per-sample calls, batch definitions)
    MAY be tracked. Each is listed in `PROJECT_CONTEXT.md` with the script that regenerates it, and is regenerated,
    never edited by hand, when its inputs change. Never for raw data, controlled-access data or personal data.
* `CODE-PROV-003` Each analysis run writes to its own directory (e.g. `results/<run>/`) containing, besides the
  outputs, a provenance file with: parameters, git commit (and whether the tree was dirty),
  package and external tool versions, input checksums or identifiers, and timestamp.
* `CODE-PROV-004` Write outputs atomically (write to a temporary file in the same directory, then rename) so a
  failed run never leaves a valid-looking partial file.
* `CODE-PROV-005` Scripts MUST NOT overwrite existing outputs unless an explicit `--overwrite` flag is given.
  Re-running a script with the same inputs and parameters SHOULD give the same outputs.
  * *Why:* re-running a script by mistake should never destroy hours of results; an explicit flag makes overwriting
    a deliberate act.

### Downloads from public archives

* `CODE-DL-001` **Every downloaded file MUST be verified against the checksum published by the archive** (e.g. ENA
  `fastq_md5` or `submitted_md5`, Synapse `contentMd5`, the `.md5` files of EGA) before it is used. This holds for
  every code path that yields a file: first download, resumed partial file, retry, parallel or multi-connection
  transfer, and a file found already present from an earlier run. A successful exit code, an HTTP 200, the expected
  size or a check performed internally by the client MUST NOT be taken as evidence of a correct file: the pipeline
  compares the checksums itself and records, per file, the observed and expected values and the verification method
  in the provenance file (`CODE-PROV-003`, run provenance).
  * If the archive publishes no checksum for a file, the download MUST fail, unless (a) the archive provides an
    integrity check for its container format, which is run and recorded instead (SRA: `prefetch --verify yes`
    followed by `vdb-validate`), or (b) the decision to accept the file without a checksum is recorded, with its
    reason, next to the source of the file (`CODE-PROV-002`, how to obtain inputs). Files derived after verification
    (e.g. `fasterq-dump` output, recompressed) record the checksum of the produced file.
  * Tables of URLs and checksums MUST be matched one to one (e.g. `zip(..., strict=True)` in Python); a list of
    checksums shorter or longer than the list of files is an error.
  * *Why:* a multi-connection client can write a file with exactly the published size and a wrong MD5; a `zip` over
    URL and checksum lists of different length drops files silently; client-side checks are skipped in some code
    paths (resume, reuse of an existing file).
* `CODE-DL-002` **Retries restart the whole unit.** MUST NOT combine a client-side retry that restarts the transfer (e.g.
  `curl --retry`) with output to a pipe (`curl … | gzip -dc | head`). Retry the whole unit (e.g. both mates of a pair)
  from the start, and validate each attempt.
  * *Why:* after a dropped connection the client sends the file again after the bytes already passed downstream; the
    consumer receives a truncated record followed by the start of the file, and a line count can still look right.
* `CODE-DL-003` **Validate streamed reads record by record**: every FASTQ record has its `@` and `+` lines and a sequence as long
  as its quality string, the expected number of records arrived, and paired streams have matching read IDs in every
  record, not only the first and the last.
  * *Note:* the first $N$ reads of a file are not a random sample of the library (they follow lane and tile order);
    state it wherever a streamed head is used as a subsample.
* `CODE-DL-004` **Treat unexpected responses as failed attempts.** Direct file URLs of archives that serve files without
  redirection (e.g. ENA over HTTPS) are fetched without following redirects (`curl` without `-L`,
  `wget --max-redirect=0`); where redirects are part of the service (e.g. signed cloud-storage URLs), check the
  content (gzip header, then checksum) before using the file.
  * *Why:* an archive backend can intermittently answer with a redirect to an HTML listing for a file that exists;
    following it stores HTML or feeds it to a decompressor.
* `CODE-DL-005` **Back off over minutes, not seconds**: several attempts per file with increasing waits (e.g. 8 attempts, waiting
  $\min(60a, 600)$ s before attempt $a$; about half an hour in total) before a batch fails.
  * *Why:* archive-side incidents can last tens of minutes, much longer than a network interruption.
* `CODE-DL-006` **Limit concurrent connections per archive and count all of them** (files in parallel times
  connections per file) against one budget per archive and source address: 3–4 by default, up to about 8 when the
  measured throughput per connection is low and the archive's terms allow it; record the budget and the measurement.
  Retry refused connections (`wget --retry-connrefused`; wget does not retry them by default).
  * *Why:* archives limit or refuse connections per source address, not per file; a multi-connection client
    multiplies the number of connections a budget counted in files suggests.
* `CODE-DL-007` **Fetch and validate the run table once** (URLs, checksums, read counts per run), track it (`CODE-PROV-002`,
  *Deviate when*), and have jobs read it instead of querying the archive API at run time.
  * *Why:* archive APIs can return empty or error responses for a fraction of accessions at any time, and
    downstream tools do not always detect them.
* `CODE-DL-008` **Test downloads from the machines that will run them.** On clusters, compute nodes can reach the internet
  through a different address than login nodes, and external services can answer differently by source.
* `CODE-DL-009` **Bound waits inside download jobs and size their time limit from measured throughput.** A job that
  waits (for a download slot, disk space or a download-ahead limit) spends its time limit while waiting: bound the
  wait, then start or fail explicitly, and add the bound to the time request. Size the transfer part from throughput
  measured on the machines that run the jobs (`CODE-DL-008`, testing from those machines;
  `HPC-JOB-002` in `HPC_GUIDELINES.md` §4, resources from measurements), add the back-off of `CODE-DL-005` (retry waits), and increase
  the time limit with the attempt number when a timed-out job is resubmitted.
  * *Why:* jobs that wait inside their allocation time out without an error of their own, some before transferring a
    byte, and the timeout looks like a slow transfer.
* `CODE-DL-010` **Choose the transfer method per repository in this order**, falling back to the next one when the
  previous is unavailable or fails repeatedly. Every method is subject to the other rules of this subsection.

  | Repository / case | Order of preference | Verification |
  |---|---|---|
  | ENA, run with FASTQ (or submitted FASTQ) and MD5 | 1. HTTPS with a multi-connection client (e.g. `aria2c -x 4 -s 4`), files of a unit in parallel within the budget of `CODE-DL-006`, whole-file retry on MD5 mismatch. 2. HTTPS with `wget` or `curl` and resume (`-c`), files in parallel. 3. FTP, as in 2 (only if HTTPS is blocked). | MD5 from the tracked run table (`CODE-DL-007`), per file |
  | SRA-only run (no FASTQ in ENA) | 1. `prefetch --verify yes` and `vdb-validate`, then `fasterq-dump --split-files --include-technical` and parallel compression. 2. The same from the cloud mirror of SRA, if reachable (test it first, `CODE-DL-008`). | `vdb-validate` on the container; checksum of the produced FASTQ recorded |
  | Synapse | 1. `synapseclient`, one client and one staging directory per file, files in parallel within the budget of `CODE-DL-006` (e.g. 3). | `contentMd5` from the tracked file table, compared by the pipeline |
  | EGA | 1. `pyega3` (or the EGA client in use), resumable. | MD5 published by EGA, from a tracked table |
  | Other HTTPS sources (GEO supplementary files, institutional servers) | 1. `wget` or `curl` with resume; a multi-connection client only if the server supports range requests. | Published checksum; if none is published, `CODE-DL-001` (a) or (b) |

  * *Why:* for large FASTQ files from ENA, multi-connection HTTPS was about 2–5 times faster than single-connection
    `wget` in a pilot (about 20–30 min versus 65–100 min per 19 GB library), and SRA adds conversion and
    recompression to the transfer. Site profiles record such measurements (e.g. `hpc_sites/DIPC_HYPERION.md`, §8).
  * *Deviate when:* the archive's terms limit connections (use the single-connection method), or measurements from
    the machines that run the jobs (`CODE-DL-008`) show another order to be faster or more reliable.


## 14. Configuration, paths and scripts

* `CODE-CONF-001` The repository MUST be runnable from a fresh clone: all configuration is tracked in git, and
  the steps from clone to results are documented in `README.md`.
  * *Why:* results that depend on untracked local files cannot be reproduced by anyone else, including the same
    person on another machine. *Trade-off:* configuration that differs between machines must be expressed through
    paths relative to the project root or environment variables.
* `CODE-CONF-002` Library code MUST NOT contain absolute or machine-specific paths. Paths enter through
  function parameters, command-line arguments or tracked configuration, and are relative to the
  project root.
* `CODE-CONF-003` In project scripts and applications, resolve the project root by locating `pyproject.toml` upwards from the
  calling file; an environment variable (e.g. `PROJECT_ROOT`) MAY override it. Library code that may be installed
  as a package does not depend on the repository layout: it receives paths as parameters, and locates its own data
  files with `importlib.resources`.
* `CODE-CONF-004` Inputs that must live outside the repository (large data) are linked into `data/` and their
  origin is documented as in §13 (data provenance).
* `CODE-CONF-005` Executable scripts put their logic in a `main()` function guarded by
  `if __name__ == "__main__":`; importing a script MUST NOT run work.
* `CODE-CONF-006` Command-line drivers use `argparse` (or `typer` if already a dependency), provide `--help`, and
  expose every analysis parameter that a user may need to change.
* `CODE-CONF-007` MUST NOT store secrets, tokens or credentials in code, configuration files or notebooks.

## 15. Reproducibility

* `CODE-REPRO-001` Expose random seeds as parameters and use local generators
  (`numpy.random.default_rng(seed)`); MUST NOT rely on implicit global random state.
* `CODE-REPRO-002` Results MUST NOT depend on the number of workers, parallel execution order, filesystem
  listing order or iteration order of unordered collections. Sort explicitly where order
  matters; derive per-task seeds from the master seed and a stable task key, so that random draws are identical.
  Floating-point results MUST agree within a stated tolerance; bit-identical results are not required when the order
  of floating-point operations legitimately changes (parallel reductions, different BLAS).
  * *Why:* parallel execution changes the order of floating-point sums and of random number draws; conclusions must
    not change whether run on a laptop or a cluster.
* `CODE-REPRO-003` Document operations whose results may vary across runs or environments.
* `CODE-REPRO-004` Keep the environment lockfile (see §3, environment) current.

## 16. Testing

* `CODE-TEST-001` Add tests for new behaviour and bug fixes. For bug fixes, add a regression test that fails
  before the fix and passes after it whenever practical.
* `CODE-TEST-002` Use pytest. Name tests after behaviour (`test_returns_empty_frame_for_empty_input`), not after
  function names.
* `CODE-TEST-003` Test edge cases and error paths, not only the happy path.
* `CODE-TEST-004` Test through public interfaces; test private helpers directly only when they contain
  substantial isolated logic. Avoid tests that mirror implementation details.
* `CODE-TEST-005` Compare numbers with tolerances (`numpy.testing.assert_allclose`) and tables with
  `polars.testing.assert_frame_equal`.
* `CODE-TEST-006` Build small input tables inline or from tiny fixtures in `tests/data/`. Tests MUST NOT read
  real project data or results, nor access the network.
  * *Why:* tests that read real data are slow, fail on machines without the data, and change their outcome when the
    data change; small fixtures make each test's expectation explicit.
* `CODE-TEST-007` Keep tests independent: no shared mutable state and no reliance on execution order.
* `CODE-TEST-008` For pipelines, SHOULD include an end-to-end test on a minimal fixture whose expected output is
  checked in.

## 17. Style

* `CODE-STYLE-001` Follow PEP 8. Formatting is handled by `ruff format`; do not hand-format against it.
* `CODE-STYLE-002` Maximum line length: 99, configured in `pyproject.toml`.
* `CODE-STYLE-003` Keep imports at module level unless a lazy or conditional import has a concrete reason.
  Group them as standard library, third-party, local (ruff rule `I`).
* `CODE-STYLE-004` Use f-strings for all string interpolation, including `logging` calls (disable ruff rule
  `G004` in `pyproject.toml`). Use `%`-style lazy logging arguments only where there is a
  measured cost, e.g. debug messages inside hot loops.
* `CODE-STYLE-005` SHOULD configure `pre-commit` with ruff (lint and format), `shellcheck` and `nbstripout`.

## 18. Dependencies

* `CODE-DEP-001` Prefer the standard library and existing project dependencies when they provide a clear
  solution; check the project before adding a library.
* `CODE-DEP-002` Add a dependency when it substantially improves correctness, maintainability, performance or
  reliability. Do not reimplement mature domain algorithms without a good reason; do not add a
  dependency for something a few clear lines can do.
* `CODE-DEP-003` Heavy or optional dependencies go into an extra in `pyproject.toml`; the core package MUST
  import without them.

## 19. Performance

* `CODE-PERF-001` Prioritize correctness and readability; optimize when there is a meaningful expected or
  measured benefit.
* `CODE-PERF-002` Avoid repeated expensive conversions, allocations, I/O and computations.
* `CODE-PERF-003` Keep optimized code understandable and document non-obvious performance decisions.
* `CODE-PERF-004` Parallelize only independent, expensive workloads where the benefit justifies the complexity
  and memory; do not parallelize small or I/O-bound work without evidence.

## 20. Logging and side effects

* `CODE-LOG-001` Use the `logging` module, not `print()`, in library code. `tqdm` progress bars are acceptable
  behind a `progress=` flag.
* `CODE-LOG-002` Functions MUST NOT write files, mutate globals or access the network unless that is their
  stated responsibility.
* `CODE-LOG-003` Avoid mutable global state. Module-level constants are encouraged; if mutable module state is
  unavoidable, make it private (`_name`), access it through functions and document why.
* `CODE-LOG-004` Prefer pure functions; isolate I/O and side effects at the edges (scripts, loaders, writers).

## 21. Compatibility and interfaces

* `CODE-IFACE-001` Treat public APIs, documented output schemas, filenames, column names and serialized formats
  as interfaces. List them in `PROJECT_CONTEXT.md` §7 (interfaces).
  * *Why:* notebooks, other repositories and people depend on these names; a rename breaks them silently or, worse,
    makes them read the wrong column.
* `CODE-IFACE-002` MUST NOT rename, remove or change the semantics of an interface unless the task requires it.
  When it must change, state the compatibility implications explicitly and record the change in
  the decisions log.
* `CODE-IFACE-003` A change to a rule that defines an output (filter, threshold, classification rule, key) makes
  previous outputs stale: say so and record it in the decisions log.
* `CODE-IFACE-004` Code shared across repositories is imported from its home package, never copied. Exception: recipes of this
  guideline repository MAY be vendored as described in §26 (recipes).

## 22. Notebooks

* `CODE-NB-001` Notebooks orchestrate and plot. Logic lives in the project package and is imported.
* `CODE-NB-002` The first code cell resolves the project root (see §14, paths) and enables
  `%load_ext autoreload` / `%autoreload 2`.
* `CODE-NB-003` Figures go to `results/figures/<notebook-prefix>_<name>.svg` (vector by default) and follow
  `FIGURE_GUIDELINES.md` when present.
* `CODE-NB-004` Notebooks MUST run top to bottom on a fresh kernel. Commit them with outputs stripped.
* `CODE-NB-005` Notebooks are excluded from ruff and exempt from the docstring rules.

## 23. Version control

* `CODE-VCS-001` Each commit is one logical change; an implementation and its tests go in the same commit.
* `CODE-VCS-002` Commit message: a summary line in the imperative, capitalized, at most 50 characters, without
  a trailing period (`Fix off-by-one in GFF parser`); then a blank line and, when useful, a body
  wrapped at 72 characters explaining why the change is needed and any side effects.
* `CODE-VCS-003` At the end of each meaningful operation or code change, propose a commit: a summary of the
  diff and the commit message.
* `CODE-VCS-004` MUST NOT commit until the user has approved the proposed diff and message. MUST NOT push
  without explicit approval.
  * *Why:* the user owns the history of the project; reviewing the diff before committing catches unintended
    changes while they are still cheap to undo.
* `CODE-VCS-005` Maintain a `.gitignore` covering data, results, environment directories and caches.
* `CODE-VCS-006` MUST NOT rewrite published history or force-push shared branches unless the user asks.

## 24. Validation

* `CODE-VAL-001` After modifying code, run the most relevant tests and static checks (`pytest`, `ruff check`,
  `ruff format --check`, `shellcheck` for shell scripts): targeted tests first, the broader
  suite when appropriate.
* `CODE-VAL-002` MUST NOT claim that tests or checks passed unless they were actually run.
* `CODE-VAL-003` If validation cannot be performed, state exactly what was not verified.

## 25. Working process (agents)

* `CODE-AGENT-001` **Ask before assuming when a request admits more than one reading that changes the result**: a
  vague instruction, a term that can name more than one thing (two files, samples, rules or versions with the same
  name), or an unstated scope. First resolve what the project answers (code, `PROJECT_CONTEXT.md`, the decisions
  log, the data); ask about what remains before doing the work that depends on it.
  * Ask one focused question per decision, with the concrete readings as options and the recommended one first;
    group independent questions in one message.
  * When one reading is clearly the most plausible, or the readings lead to the same result, proceed and state the
    assumption explicitly in the reply, next to the result it affects, not only in code or logs. Assumptions that
    outlive the task go to the decisions log.
  * Work that cannot be undone (`CODE-AGENT-004`, deletion and overwriting; `CODE-AGENT-005`, long or expensive
    jobs; pushes to shared branches) is never done under an unconfirmed assumption.
  * *Why:* an agent that silently picks one reading delivers a coherent answer to a different question, and the
    error surfaces only when someone checks the premise; a stated assumption costs one sentence and can be corrected
    at once.
  * *Trade-off:* questions interrupt the user; asking what the project already answers, or about choices that do not
    change the result, moves the agent's work to the user.
  * *Deviate when:* no one can answer in time (unattended runs): choose the most plausible reversible reading, record
    it as an assumption in the output, and flag it for review.
* `CODE-AGENT-002` State important assumptions explicitly rather than silently choosing defaults, especially
  analysis parameters, thresholds and reference versions (`CODE-AGENT-001`, stating assumptions).
* `CODE-AGENT-003` MUST NOT fabricate results, numbers, citations, accessions or tool outputs. Report only what
  was computed or read, and say where it came from.
* `CODE-AGENT-004` MUST NOT delete or overwrite data, results or user files without explicit confirmation.
* `CODE-AGENT-005` Ask before launching jobs that are long-running, expensive, or download large resources.
  * *Why:* compute, storage and shared resources are limited and their use is the user's decision; a misconfigured
    long job can also occupy resources for hours before the error is visible.
* `CODE-AGENT-006` Do not create documents, summaries or example scripts beyond the documentation set of §2
  (documentation set) unless asked.
* `CODE-AGENT-007` When presenting edits to an existing file, show only the changed sections unless the file is
  short.

## 26. Related guidelines

Available in the same repository; copy the ones relevant to the project:

* `STATISTICS_GUIDELINES.md` — experimental design, test choice, multiple testing, reporting.
* `FIGURE_GUIDELINES.md` — chart choice, colour, layout, export (planned; not yet available).
* `RNASEQ_BULK_GUIDELINES.md` — bulk RNA-seq from reads to differential expression and enrichment.
* `SCRNASEQ_GUIDELINES.md` — single-cell and single-nucleus RNA-seq.
* `SPATIAL_GUIDELINES.md` — spatial transcriptomics (Visium, Xenium, MERFISH: multiplexed error-robust fluorescence in situ hybridization).
* `META_ANALYSIS_GUIDELINES.md` — combining results across studies, including omics data.
* `SEX_INFERENCE_GUIDELINES.md` — sex-chromosome complement from expression; reference routing.
* `SAMPLE_IDENTITY_GUIDELINES.md` — duplicates, swaps, mixtures and genotype fingerprints.
* `DECONVOLUTION_GUIDELINES.md` — cell-type composition of bulk samples.
* `HPC_GUIDELINES.md` — shared clusters with a batch scheduler: storage roles, data flow, containers, jobs; site
  facts in `hpc_sites/`.
* `STATISTICS_REFERENCES.md`, `OMICS_REFERENCES.md` — references cited by the guidelines above (`references/` upstream,
  `guidelines/references/` in a project).

### Recipes (reference implementations)

* `CODE-RCP-001` Non-obvious procedures of the domain guidelines MAY have a reference implementation in `recipes/` of the same
  repository, linked from the guideline as `Recipe: recipes/<area>/<module>.py`. `recipes/README.md` lists, per
  guideline area, which procedures are implemented, delegated to an external tool, or planned, with the scope of
  their tests. Each recipe states the guideline section it implements, its inputs and outputs, and its defaults
  with their source; tests on synthetic data are in `recipes/tests/`. Tests on synthetic data check the contract of
  the code, not the validity of the method on real data.
* `CODE-RCP-002` **When a guideline links a recipe, use it instead of re-implementing the procedure**. Preferred: install the
  repository at a pinned commit and import it. Alternative (vendoring): copy the module together with every
  `recipes` module it imports, its tests, the license and the source commit. Record which one, and the commit, in
  `PROJECT_CONTEXT.md`.
  * *Why:* several agents implementing the same procedure from prose diverge in details (seeds, subsampling, how
    clusters are matched, which tail of a distribution is filtered), and their results stop being comparable.
* `CODE-RCP-003` Adapt a recipe through its parameters and named constants. A change to its logic is a deviation: record it under
  *Overrides* in `PROJECT_CONTEXT.md` and keep the recipe's tests passing (or replace them with tests of the new
  behaviour).
* `CODE-RCP-004` If a recipe and its guideline disagree, the guideline wins; report the discrepancy upstream.
