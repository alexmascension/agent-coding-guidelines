# Coding guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 2.3.0 (2026-10-14)
>
> The key words **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT** and **MAY** are used as in
> RFC 2119. MUST rules are hard constraints. SHOULD rules are strong defaults; deviating from
> them requires a reason stated in the change or in `PROJECT_CONTEXT.md`.

## 0. How to use this document

* This file is project-agnostic. Copy it **verbatim** into the root of every project.
* Everything specific to a project (stack details, layout, vocabulary, data invariants,
  conventions, interfaces, decisions) goes into `PROJECT_CONTEXT.md` (see §2, documentation
  set).
* Do not edit the copy of this file to adapt it to a project. Record adaptations in the
  *Overrides* section of `PROJECT_CONTEXT.md`. This keeps the copy comparable with upstream and
  lets it be updated without losing project decisions.
* When this file and `PROJECT_CONTEXT.md` conflict, `PROJECT_CONTEXT.md` wins **only** for
  rules it explicitly lists under *Overrides*. Otherwise this file wins.
* Domain guidelines from the same repository (statistics, figures, RNA-seq, …) are copied
  alongside this file when relevant to the project and follow the same override mechanism.
* **Versioning.** Every guideline file carries a semantic version, and changes are recorded in
  the repository `CHANGELOG.md`:
  * MAJOR: a section is renumbered or removed, or a rule is removed or reversed (existing
    overrides may point to the wrong rule).
  * MINOR: new rules or sections appended without renumbering.
  * PATCH: wording or clarifications that do not change what a rule requires.

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
| `CODING_GUIDELINES.md` | This file, unmodified. |
| `README.md` | What the project does, how to install it, how to run it, where outputs go. |

Add further documents only when they carry information that does not fit the above, e.g.
`PROJECT_SCOPE.md` (scientific scope and design rationale), `DECISIONS.md` (if the decisions log
outgrows `PROJECT_CONTEXT.md`), copied domain guidelines (`STATISTICS_GUIDELINES.md`,
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
2. `CODING_GUIDELINES.md` — general development rules.
3. Domain guidelines present in the root (`*_GUIDELINES.md`), when the task touches them.
```

For agents that do not read `AGENTS.md` (e.g. Claude Code reads `CLAUDE.md`), add the file the
agent expects containing only a reference to `AGENTS.md` (in `CLAUDE.md`: `@AGENTS.md`).

### Required structure of `PROJECT_CONTEXT.md`

`PROJECT_CONTEXT.md` MUST keep the following sections, in this order. A section that does not
apply MUST be kept with the text `None.` so its absence is an explicit decision rather than an
omission.

1. **Purpose** — one paragraph; link to `README.md` / `PROJECT_SCOPE.md` for detail.
2. **Stack** — languages, key libraries, external tools **with versions**, environment name and
   how to create it.
3. **Layout** — directories, what goes in each, where inputs and outputs live, which paths
   these guidelines apply to and which are exempt.
4. **Domain vocabulary** — controlled labels, categorical values and column names, and the
   module where they are declared (see §6, controlled vocabulary).
5. **Data model and invariants** — atomic tables, entity keys, derivation rules, filters that
   define each table (see §7, tables and schemas).
6. **Conventions** — coordinate systems, genome build and annotation release, units,
   normalization, thresholds and their rationale (see §11, numerical code).
7. **Interfaces** — outputs consumed by notebooks, other repositories or people; sibling
   libraries that MUST be reused rather than copied (see §21, interfaces).
8. **Overrides** — first, the version of each adopted guideline file; then every deviation
   from them, each with the rule reference, the replacement rule and the reason.
9. **Decisions log** — dated entries: decision, rationale, affected outputs. MAY instead link to
   `DECISIONS.md`.

Keep `PROJECT_CONTEXT.md` synchronized with the code. A change that alters a convention, an
interface, a threshold or a data invariant MUST update it in the same change.

## 3. Languages and tooling

* Python 3.11+ is the default language.
* Prefer a maintained Python implementation of an established method over calling R or another
  language. Use it only if it is validated against the reference implementation; when results
  must match the reference (published method, comparison with prior work), state which
  implementation and version was used, and document known differences.
  * *Why:* one language keeps the environment, the data structures and the tests in one place, and Python code is
    easier for the project's readers to inspect. *Trade-off:* ports can differ numerically from the reference;
    hence the validation requirement. *Deviate when:* no port exists, or the port is not validated for the analysis
    at hand.
* Use R, bash or other languages only when no adequate Python option exists or when the user
  requests it. Isolate such code in standalone scripts with a file-based interface (inputs and
  outputs on disk, parameters on the command line) called from Python. Avoid in-process bridges
  such as `rpy2` unless there is a concrete reason.
  * *Why:* a file interface makes every intermediate result inspectable and lets each language run in its own
    environment; in-process bridges tie the versions of two runtimes together and fail in ways that are hard to
    reproduce.
* Bash is for glue only (calling tools, moving files, little data manipulation):
  * Rewrite a script in Python once it exceeds about 100 lines or needs non-trivial control
    flow.
  * Start every script with `set -euo pipefail`; quote all variable expansions.
  * Use `$(…)` rather than backticks and `[[ … ]]` rather than `[ … ]`; send errors to STDERR;
    put the logic in a `main` function called at the end of the file.
  * Scripts MUST pass `shellcheck`.
* R code, when needed, follows the tidyverse style guide, is formatted with `styler` and linted
  with `lintr`.
* Default Python stack (a project MAY override it in `PROJECT_CONTEXT.md`): polars for tables,
  NumPy/SciPy for numerics, `statsmodels` for statistical models, matplotlib/seaborn for figures,
  pytest for tests, ruff for linting and formatting.
* Environment: a conda environment plus `pip install -e .` of the project package. Export the
  resolved environment (`conda env export --no-builds > environment.lock.yml` or `conda-lock`)
  whenever dependencies change.
  * *Why:* an environment file with version ranges resolves to different packages months later; the lockfile
    records what actually ran.
* Default layout (a project MAY override it):
  `src/<package>/` library code, `src/scripts/` command-line drivers, `src/notebooks/`
  analysis notebooks, `tests/`, `data/` inputs, `results/` outputs.
  These guidelines apply to library code, scripts and tests; notebooks follow §22 (notebooks).

## 4. General principles

* Write clean, readable and maintainable code.
* Prefer simple, explicit solutions over clever or unnecessarily complex ones.
* Give each function one clear responsibility; split large or deeply nested logic into focused
  helpers when this improves clarity. A function longer than about 40 lines is a signal to
  consider splitting it, not a hard limit.
* Extract repeated logic into reusable helpers; do not introduce abstractions unless they
  genuinely simplify the code.
* Use descriptive names and intermediate variables when they make logic easier to follow.
  Avoid compact one-liners that reduce readability; a comprehension SHOULD have at most one
  `for` clause and one condition.

## 5. Changes to existing code

* Preserve existing behaviour and interfaces unless the requested change requires otherwise.
* Make the smallest coherent change that fully solves the requested problem.
  * *Why:* small diffs can be reviewed; unrelated changes mixed into a fix hide regressions and make reverting one
    change impossible without reverting the other.
* MUST NOT refactor or rewrite unrelated code while implementing a feature or bug fix.
* Reuse existing project utilities and conventions instead of creating parallel implementations.
* If a broader refactor would materially improve the solution, propose it separately.
* If you notice an unrelated bug, report it; do not fix it in the same change.

## 6. Naming and controlled vocabulary

* `UPPER_SNAKE_CASE` for module-level and configuration constants; `lowercase_with_underscores`
  for functions, parameters, variables and other non-constant names.
* Prefer descriptive, self-contained names over abbreviations, unless the abbreviation is
  standard in the domain (`tpm`, `fdr`, `cds`).
* Loop and comprehension variables stay descriptive when their meaning is not obvious; avoid
  `x`, `a`, `tmp`.
* Short names are allowed only for standard mathematical notation (`n`, `k`, `p`, `X`), for
  notation that follows a reference paper or algorithm (cite the source in the docstring), or
  for labels defined as such in `PROJECT_CONTEXT.md`. Public APIs keep descriptive names.
* For hierarchical information, order name components from broader to more specific
  (`gene_copy_count`, not `count_copy_gene`).
* **Controlled vocabulary has a single source of truth.** Categorical labels, region names,
  classification values and column names MUST be declared once (a `config.py` or `schema.py`
  module, as listed in `PROJECT_CONTEXT.md` §4, domain vocabulary) and imported everywhere
  else. MUST NOT invent synonyms, aliases or alternative spellings in modules or notebooks.
  * *Why:* two spellings of one label (`XAR` and `X_added`) silently split one group into two, and every downstream
    count is wrong without any error being raised.

## 7. Data, tables and schemas

* polars is the default for every table. pandas appears only at boundaries that require it
  (plotting libraries, `statsmodels`, third-party APIs): convert once, as late as possible, and
  never round-trip.
* MUST NOT migrate existing code between pandas and polars solely for stylistic reasons.
* Schemas are part of the function contract: document required columns and important dtypes
  (see §9, documentation), and validate them at public entry points.
* Preserve row order, schemas and dtypes unless the operation explicitly changes them. In
  polars, `group_by` does not preserve group order; use `maintain_order=True` or sort explicitly
  when order matters downstream.
* Joins MUST state their expected cardinality (polars `validate="1:1"` / `"m:1"`) when a
  violation would silently duplicate or drop rows, and SHOULD check row counts before and after.
  * *Why:* a many-to-many join on a key that was assumed unique duplicates rows without error, inflating every
    count computed afterwards.
* Prefer columnar or vectorized operations for bulk transformations; do not force
  vectorization when row-oriented logic is clearer or required by an external API.
* **Aggregates are derived, never accumulated.** Keep the most granular (atomic) table as the
  source of truth. Re-derive summaries from it after any filter; never obtain them by combining
  previously aggregated rows.
  * *Why:* summing pre-aggregated rows double counts entities that appear in several groups and keeps results
    computed under old filters; re-deriving from the atomic table guarantees that every summary reflects the
    current filters.
* **Entity keys are stable.** Each entity has a documented key (`PROJECT_CONTEXT.md` §5, data
  model) that is identical across modules. Any rule that merges or splits entities is a named
  constant and is applied in one place.
* Use Parquet for intermediate and output tables (it preserves dtypes). Use TSV/CSV only for
  human-facing exports or tool interoperability.
  * *Why:* CSV loses dtypes: identifiers with leading zeros become numbers, nulls become empty strings, and
    categorical columns are re-parsed on every load.

## 8. Type hints

* Add type hints to public functions and non-trivial internal functions.
* Use built-in generics (`list[str]`, `dict[str, int]`) and `X | None`, not `typing.List` or
  `Optional[X]`.
* Do not add annotations that reduce readability or that cannot express the type meaningfully.

## 9. Documentation

### Language and docstrings

* Write repository documentation, notebook prose, code comments, docstrings, log and error
  messages, and figure/table labels in **English**, regardless of the conversation language.
* Write NumPy-style docstrings for public functions and non-trivial internal functions:

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

* For functions that take or return tables, the docstring MUST list the required columns and
  their dtypes, and state any ordering or uniqueness assumptions.
* Document assumptions, invariants, edge cases and non-obvious behaviour.

### Self-contained comments and documents

* Every comment and document MUST be understandable by a reader who knows the project scope but
  has no knowledge of how the code was developed.
* Comments explain *why*, not *what*. Use them for non-obvious decisions, assumptions,
  performance considerations and tricky logic; do not comment self-explanatory code.
* State the **reason**, not the **history**. MUST NOT refer to the conversation, request,
  meeting, previous version or decision that led to the code (`as discussed`, `per user
  request`, `changed after review`, `previously we used X`). Write the reason itself
  (`# Merge hits closer than MERGE_HITS_BP: aligners split one copy at short indels`).
  Decision history belongs only in the decisions log.
* MUST NOT leave commented-out code; version control keeps old code.
* `TODO` comments use the format `# TODO: <issue or link> - <description>` and are reported to
  the user when created.
* Documents that define rules (guidelines, `PROJECT_CONTEXT.md`) give each non-obvious rule its rationale: the
  mechanism (*Why*), its cost (*Trade-off*) and the justified exceptions (*Deviate when*).
  * *Why:* a rule whose reason is visible can be checked and applied to cases it does not mention; a bare rule
    can only be obeyed or broken.
* Keep documentation minimal and current: delete or update documents that no longer describe
  the code. Outdated documentation is worse than missing documentation.
* Keep documentation synchronized with the implementation, including `README.md` and
  `PROJECT_CONTEXT.md` when an interface or design decision changes.

### Markdown

* One H1 per file; ATX headings (`#`), each heading unique within the file.
* Fenced code blocks always declare their language (`python`, `bash`, `text`).
* Inline equations use `$...$`; display equations use `$$...$$` on their own lines.
  MUST NOT use `\(...\)` or `\[...\]`. This also applies to notebook Markdown cells.
* **Cross-references state what they point to.** A reference to a section, rule, figure, table
  or numbered item MUST be followed by a short description, e.g. `§13 (data provenance)`,
  `Fig. 2 (copy-age distribution)`, `step 4 (normalization)`. The description need not be the
  exact title, only enough to identify the target without opening it.

## 10. Robustness

* Validate important inputs and fail early with error messages that precisely match the failing
  condition, show interpolated values unambiguously (`f"Not a probability: {p=}"`) and are easy
  to grep.
* MUST NOT use `assert` to validate inputs or enforce program logic: assertions are removed
  under `python -O`. Use `assert` only in tests and for internal invariants whose removal would
  not change behaviour.
* Handle expected edge cases explicitly (empty inputs, missing groups, single-sample designs).
* MUST NOT silently swallow errors. Avoid broad `except Exception` unless there is a specific,
  documented reason, and re-raise or log with context.
* MUST NOT use mutable objects (`[]`, `{}`, a DataFrame) as default argument values; use `None`
  and create the object inside the function.
  * *Why:* Python evaluates default values once, at function definition, so a mutable default is shared by all
    calls and accumulates state between them.
* MUST NOT mutate input arguments (DataFrames, arrays, lists, dicts) in place unless this is
  intentional and documented.
* Open files and other resources with `with` statements.
* Avoid "power features" (metaclasses, monkeypatching, dynamic attribute access via
  `getattr`/`setattr` on computed names, import-time side effects) unless there is no simpler
  alternative; document why when used.

## 11. Numerical and scientific code

* Make units, coordinate systems, normalization conventions and statistical assumptions explicit
  when they are not obvious.
* Avoid exact floating-point equality unless mathematically guaranteed; use tolerances.
* Handle NaN, null, infinity and empty inputs deliberately. In polars, `null` and `NaN` are
  distinct; decide which one each column uses.
* Check for `None` with `is None` / `is not None`, never by truthiness: `0`, `0.0` and empty
  containers are falsy, and the truth value of an array or DataFrame is ambiguous or raises.
* MUST NOT silently change numerical precision or relevant dtypes.
* Prefer established numerical, statistical and domain implementations over custom ones.
* **One internal coordinate convention.** Use a single internal convention (default: 0-based,
  half-open). Convert every external format at the parser, at load time, never later. Record the
  convention of each input format in `PROJECT_CONTEXT.md` §6 (conventions).
  * *Why:* off-by-one errors arise at every conversion between 0- and 1-based systems; converting once, at the
    parser, confines the risk to one tested function.
* Genome build and annotation release are part of every coordinate-bearing table and are
  checked when tables are combined.
* **Every reported number states its universe**: which input set, which method or branch, and
  which filters it was computed on. A count without its universe is not a result.
* **Classifications carry their rule.** A column produced by a documented classification rule
  (as opposed to a direct observation) MUST be reported together with the identifier or version
  of the rule that produced it.
  * *Why:* the same column name produced by two versions of a rule looks comparable but is not; without the rule
    version, results from different runs cannot be reconciled.
* Thresholds and analysis parameters are named constants in one configuration module, each with
  a short comment giving its rationale or source. MUST NOT hard-code them inline.
  * *Why:* inline values get copied and then diverge between modules; a named constant with its rationale is
    changed in one place and is visible in review.
* Statistical practice (test choice, multiple testing, effect sizes, reporting) follows
  `STATISTICS_GUIDELINES.md` when present.

## 12. External tools and subprocesses

* Call external tools with `subprocess.run([...], check=True)` using an argument list.
  MUST NOT use `shell=True` with interpolated values.
  * *Why:* with `shell=True`, file names with spaces or special characters are split or interpreted by the shell,
    and interpolated values can execute arbitrary commands.
* Check that required executables exist (`shutil.which`) before starting long work, and fail
  with a message naming the missing tool.
* Log the exact command line and the tool version for every external call that produces a
  result; store the version with the outputs (see §13, provenance).
* Capture stderr and include its tail in the error message when a tool fails.
* Parse tool outputs in dedicated parser functions that convert to internal conventions (see
  §11, coordinate convention) and return typed tables.

## 13. Data management and provenance

* Raw inputs (`data/raw/` or equivalent) are read-only. MUST NOT modify, overwrite or delete
  them.
* MUST NOT commit data, results or large binaries to git. Record how to obtain each input
  (URL, accession, release, checksum) in `README.md` or a download script.
* Each analysis run writes to its own directory (e.g. `results/<run>/`) containing, besides the
  outputs, a provenance file with: parameters, git commit (and whether the tree was dirty),
  package and external tool versions, input checksums or identifiers, and timestamp.
* Write outputs atomically (write to a temporary file in the same directory, then rename) so a
  failed run never leaves a valid-looking partial file.
* Scripts MUST NOT overwrite existing outputs unless an explicit `--overwrite` flag is given.
  Re-running a script with the same inputs and parameters SHOULD give the same outputs.
  * *Why:* re-running a script by mistake should never destroy hours of results; an explicit flag makes overwriting
    a deliberate act.

## 14. Configuration, paths and scripts

* The repository MUST be runnable from a fresh clone: all configuration is tracked in git, and
  the steps from clone to results are documented in `README.md`.
  * *Why:* results that depend on untracked local files cannot be reproduced by anyone else, including the same
    person on another machine. *Trade-off:* configuration that differs between machines must be expressed through
    paths relative to the project root or environment variables.
* Library code MUST NOT contain absolute or machine-specific paths. Paths enter through
  function parameters, command-line arguments or tracked configuration, and are relative to the
  project root.
* Resolve the project root by locating `pyproject.toml` upwards from the calling file; an
  environment variable (e.g. `PROJECT_ROOT`) MAY override it.
* Inputs that must live outside the repository (large data) are linked into `data/` and their
  origin is documented as in §13 (data provenance).
* Executable scripts put their logic in a `main()` function guarded by
  `if __name__ == "__main__":`; importing a script MUST NOT run work.
* Command-line drivers use `argparse` (or `typer` if already a dependency), provide `--help`, and
  expose every analysis parameter that a user may need to change.
* MUST NOT store secrets, tokens or credentials in code, configuration files or notebooks.

## 15. Reproducibility

* Expose random seeds as parameters and use local generators
  (`numpy.random.default_rng(seed)`); MUST NOT rely on implicit global random state.
* Results MUST NOT depend on the number of workers, parallel execution order, filesystem
  listing order or iteration order of unordered collections. Sort explicitly where order
  matters; derive per-task seeds from the master seed and a stable task key.
  * *Why:* parallel execution changes the order of floating-point sums and of random number draws; results must be
    identical whether run on a laptop or a cluster.
* Document operations whose results may vary across runs or environments.
* Keep the environment lockfile (see §3, environment) current.

## 16. Testing

* Add tests for new behaviour and bug fixes. For bug fixes, add a regression test that fails
  before the fix and passes after it whenever practical.
* Use pytest. Name tests after behaviour (`test_returns_empty_frame_for_empty_input`), not after
  function names.
* Test edge cases and error paths, not only the happy path.
* Test through public interfaces; test private helpers directly only when they contain
  substantial isolated logic. Avoid tests that mirror implementation details.
* Compare numbers with tolerances (`numpy.testing.assert_allclose`) and tables with
  `polars.testing.assert_frame_equal`.
* Build small input tables inline or from tiny fixtures in `tests/data/`. Tests MUST NOT read
  real project data or results, nor access the network.
  * *Why:* tests that read real data are slow, fail on machines without the data, and change their outcome when the
    data change; small fixtures make each test's expectation explicit.
* Keep tests independent: no shared mutable state and no reliance on execution order.
* For pipelines, SHOULD include an end-to-end test on a minimal fixture whose expected output is
  checked in.

## 17. Style

* Follow PEP 8. Formatting is handled by `ruff format`; do not hand-format against it.
* Maximum line length: 99, configured in `pyproject.toml`.
* Keep imports at module level unless a lazy or conditional import has a concrete reason.
  Group them as standard library, third-party, local (ruff rule `I`).
* Use f-strings for all string interpolation, including `logging` calls (disable ruff rule
  `G004` in `pyproject.toml`). Use `%`-style lazy logging arguments only where there is a
  measured cost, e.g. debug messages inside hot loops.
* SHOULD configure `pre-commit` with ruff (lint and format), `shellcheck` and `nbstripout`.

## 18. Dependencies

* Prefer the standard library and existing project dependencies when they provide a clear
  solution; check the project before adding a library.
* Add a dependency when it substantially improves correctness, maintainability, performance or
  reliability. Do not reimplement mature domain algorithms without a good reason; do not add a
  dependency for something a few clear lines can do.
* Heavy or optional dependencies go into an extra in `pyproject.toml`; the core package MUST
  import without them.

## 19. Performance

* Prioritize correctness and readability; optimize when there is a meaningful expected or
  measured benefit.
* Avoid repeated expensive conversions, allocations, I/O and computations.
* Keep optimized code understandable and document non-obvious performance decisions.
* Parallelize only independent, expensive workloads where the benefit justifies the complexity
  and memory; do not parallelize small or I/O-bound work without evidence.

## 20. Logging and side effects

* Use the `logging` module, not `print()`, in library code. `tqdm` progress bars are acceptable
  behind a `progress=` flag.
* Functions MUST NOT write files, mutate globals or access the network unless that is their
  stated responsibility.
* Avoid mutable global state. Module-level constants are encouraged; if mutable module state is
  unavoidable, make it private (`_name`), access it through functions and document why.
* Prefer pure functions; isolate I/O and side effects at the edges (scripts, loaders, writers).

## 21. Compatibility and interfaces

* Treat public APIs, documented output schemas, filenames, column names and serialized formats
  as interfaces. List them in `PROJECT_CONTEXT.md` §7 (interfaces).
  * *Why:* notebooks, other repositories and people depend on these names; a rename breaks them silently or, worse,
    makes them read the wrong column.
* MUST NOT rename, remove or change the semantics of an interface unless the task requires it.
  When it must change, state the compatibility implications explicitly and record the change in
  the decisions log.
* A change to a rule that defines an output (filter, threshold, classification rule, key) makes
  previous outputs stale: say so and record it in the decisions log.
* Code shared across repositories is imported from its home package, never copied.

## 22. Notebooks

* Notebooks orchestrate and plot. Logic lives in the project package and is imported.
* The first code cell resolves the project root (see §14, paths) and enables
  `%load_ext autoreload` / `%autoreload 2`.
* Figures go to `results/figures/<notebook-prefix>_<name>.svg` (vector by default) and follow
  `FIGURE_GUIDELINES.md` when present.
* Notebooks MUST run top to bottom on a fresh kernel. Commit them with outputs stripped.
* Notebooks are excluded from ruff and exempt from the docstring rules.

## 23. Version control

* Each commit is one logical change; an implementation and its tests go in the same commit.
* Commit message: a summary line in the imperative, capitalized, at most 50 characters, without
  a trailing period (`Fix off-by-one in GFF parser`); then a blank line and, when useful, a body
  wrapped at 72 characters explaining why the change is needed and any side effects.
* At the end of each meaningful operation or code change, propose a commit: a summary of the
  diff and the commit message.
* MUST NOT commit until the user has approved the proposed diff and message. MUST NOT push
  without explicit approval.
  * *Why:* the user owns the history of the project; reviewing the diff before committing catches unintended
    changes while they are still cheap to undo.
* Maintain a `.gitignore` covering data, results, environment directories and caches.
* MUST NOT rewrite published history or force-push shared branches unless the user asks.

## 24. Validation

* After modifying code, run the most relevant tests and static checks (`pytest`, `ruff check`,
  `ruff format --check`, `shellcheck` for shell scripts): targeted tests first, the broader
  suite when appropriate.
* MUST NOT claim that tests or checks passed unless they were actually run.
* If validation cannot be performed, state exactly what was not verified.

## 25. Working process (agents)

* If ambiguity materially affects correctness, interfaces, data or expected behaviour, ask
  before implementing. For minor ambiguities, choose the simplest reasonable interpretation and
  state the assumption.
* State important assumptions explicitly rather than silently choosing defaults, especially
  analysis parameters, thresholds and reference versions.
* MUST NOT fabricate results, numbers, citations, accessions or tool outputs. Report only what
  was computed or read, and say where it came from.
* MUST NOT delete or overwrite data, results or user files without explicit confirmation.
* Ask before launching jobs that are long-running, expensive, or download large resources.
  * *Why:* compute, storage and shared resources are limited and their use is the user's decision; a misconfigured
    long job can also occupy resources for hours before the error is visible.
* Do not create documents, summaries or example scripts beyond the documentation set of §2
  (documentation set) unless asked.
* When presenting edits to an existing file, show only the changed sections unless the file is
  short.

## 26. Related guidelines

Available in the same repository; copy the ones relevant to the project:

* `STATISTICS_GUIDELINES.md` — experimental design, test choice, multiple testing, reporting.
* `FIGURE_GUIDELINES.md` — chart choice, colour, layout, export.
* `RNASEQ_BULK_GUIDELINES.md` — bulk RNA-seq from reads to differential expression and enrichment.
* `SCRNASEQ_GUIDELINES.md` — single-cell and single-nucleus RNA-seq.
* `SPATIAL_GUIDELINES.md` — spatial transcriptomics (Visium, Xenium, MERFISH: multiplexed error-robust fluorescence in situ hybridization).
* `META_ANALYSIS_GUIDELINES.md` — combining results across studies, including omics data.
* `SEX_INFERENCE_GUIDELINES.md` — sex-chromosome complement from expression; reference routing.
* `SAMPLE_IDENTITY_GUIDELINES.md` — duplicates, swaps, mixtures and genotype fingerprints.
* `DECONVOLUTION_GUIDELINES.md` — cell-type composition of bulk samples.
* `STATISTICS_REFERENCES.md`, `OMICS_REFERENCES.md` — references cited by the guidelines above.

### Recipes (reference implementations)

* Non-obvious procedures of the domain guidelines have a tested reference implementation in `recipes/` of the same
  repository, linked from the guideline as `Recipe: recipes/<area>/<module>.py`. Each recipe states the guideline
  section it implements, its inputs and outputs, and its defaults with their source; tests on synthetic data are in
  `recipes/tests/`.
* **When a guideline links a recipe, use it instead of re-implementing the procedure**: copy the module (with its
  tests) into the project, or install the repository at a pinned commit and import it. Record which one, and the
  commit, in `PROJECT_CONTEXT.md`.
  * *Why:* several agents implementing the same procedure from prose diverge in details (seeds, subsampling, how
    clusters are matched, which tail of a distribution is filtered), and their results stop being comparable.
* Adapt a recipe through its parameters and named constants. A change to its logic is a deviation: record it under
  *Overrides* in `PROJECT_CONTEXT.md` and keep the recipe's tests passing (or replace them with tests of the new
  behaviour).
* If a recipe and its guideline disagree, the guideline wins; report the discrepancy upstream.
