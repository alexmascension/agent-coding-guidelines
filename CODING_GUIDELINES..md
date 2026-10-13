# Coding guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 2.0.0 (2026-10-13)
>
> The key words **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT** and **MAY** are used as in
> RFC 2119. MUST rules are hard constraints. SHOULD rules are strong defaults; deviating from
> them requires a reason stated in the change or in `PROJECT_CONTEXT.md`.

## 0. How to use this document

* This file is project-agnostic. Copy it **verbatim** into the root of every project.
* Everything specific to a project (stack details, layout, vocabulary, data invariants,
  conventions, interfaces, decisions) goes into `PROJECT_CONTEXT.md` (see §2).
* Do not edit the copy of this file to adapt it to a project. Record adaptations in the
  *Overrides* section of `PROJECT_CONTEXT.md`. This keeps the copy comparable with upstream and
  lets it be updated without losing project decisions.
* When this file and `PROJECT_CONTEXT.md` conflict, `PROJECT_CONTEXT.md` wins **only** for
  rules it explicitly lists under *Overrides*. Otherwise this file wins.
* Domain guidelines from the same repository (statistics, figures, RNA-seq, …) are copied
  alongside this file when relevant to the project and follow the same override mechanism.

## 1. Priority order

When rules or goals conflict, prioritize:

**correctness > readability > backward compatibility > consistency with the existing codebase >
brevity > performance**

## 2. Project documentation set

Every project developed with agents MUST contain at least:

| File | Content |
|---|---|
| `PROJECT_CONTEXT.md` | Project-specific decisions and constraints for development (structure below). |
| `CODING_GUIDELINES.md` | This file, unmodified. |
| `README.md` | What the project does, how to install it, how to run it, where outputs go. |

Add further documents only when they carry information that does not fit the above, e.g.
`PROJECT_SCOPE.md` (scientific scope and design rationale), `DECISIONS.md` (if the decisions log
outgrows `PROJECT_CONTEXT.md`), or copied domain guidelines (`STATISTICS_GUIDELINES.md`,
`FIGURE_GUIDELINES.md`, `RNASEQ_GUIDELINES.md`, …).

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
   module where they are declared (see §6).
5. **Data model and invariants** — atomic tables, entity keys, derivation rules, filters that
   define each table (see §7).
6. **Conventions** — coordinate systems, genome build and annotation release, units,
   normalization, thresholds and their rationale (see §11).
7. **Interfaces** — outputs consumed by notebooks, other repositories or people; sibling
   libraries that MUST be reused rather than copied (see §21).
8. **Overrides** — every deviation from this file or from a copied domain guideline, each with
   the rule number, the replacement rule and the reason.
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
* Use R, bash or other languages only when no adequate Python option exists or when the user
  requests it. Isolate such code in standalone scripts with a file-based interface (inputs and
  outputs on disk, parameters on the command line) called from Python. Avoid in-process bridges
  such as `rpy2` unless there is a concrete reason.
* Bash is for glue only. Bash scripts MUST start with `set -euo pipefail` and quote all variable
  expansions.
* Default Python stack (a project MAY override it in `PROJECT_CONTEXT.md`): polars for tables,
  NumPy/SciPy for numerics, `statsmodels` for statistical models, matplotlib/seaborn for figures,
  pytest for tests, ruff for linting and formatting.
* Environment: a conda environment plus `pip install -e .` of the project package. Export the
  resolved environment (`conda env export --no-builds > environment.lock.yml` or `conda-lock`)
  whenever dependencies change.
* Default layout (a project MAY override it):
  `src/<package>/` library code, `src/scripts/` command-line drivers, `src/notebooks/`
  analysis notebooks, `tests/`, `data/` inputs, `results/` outputs.
  These guidelines apply to library code, scripts and tests; notebooks follow §22.

## 4. General principles

* Write clean, readable and maintainable code.
* Prefer simple, explicit solutions over clever or unnecessarily complex ones.
* Give each function one clear responsibility; split large or deeply nested logic into focused
  helpers when this improves clarity.
* Extract repeated logic into reusable helpers; do not introduce abstractions unless they
  genuinely simplify the code.
* Use descriptive names and intermediate variables when they make logic easier to follow.
  Avoid compact one-liners that reduce readability.
* Preserve existing behaviour and interfaces unless the requested change requires otherwise.

## 5. Changes to existing code

* Make the smallest coherent change that fully solves the requested problem.
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
  `x`, `a`, `tmp`. Single-letter names are allowed only for standard mathematical or domain
  notation (`n`, `k`, `p`, `X`) or for labels defined as such in `PROJECT_CONTEXT.md`.
* For hierarchical information, order name components from broader to more specific
  (`gene_copy_count`, not `count_copy_gene`).
* **Controlled vocabulary has a single source of truth.** Categorical labels, region names,
  classification values and column names MUST be declared once (a `config.py` or `schema.py`
  module, as listed in `PROJECT_CONTEXT.md` §4) and imported everywhere else. MUST NOT invent
  synonyms, aliases or alternative spellings in modules or notebooks.

## 7. Data, tables and schemas

* polars is the default for every table. pandas appears only at boundaries that require it
  (plotting libraries, `statsmodels`, third-party APIs): convert once, as late as possible, and
  never round-trip.
* MUST NOT migrate existing code between pandas and polars solely for stylistic reasons.
* Schemas are part of the function contract: document required columns and important dtypes
  (§9), and validate them at public entry points.
* Preserve row order, schemas and dtypes unless the operation explicitly changes them. In
  polars, `group_by` does not preserve group order; use `maintain_order=True` or sort explicitly
  when order matters downstream.
* Joins MUST state their expected cardinality (polars `validate="1:1"` / `"m:1"`) when a
  violation would silently duplicate or drop rows, and SHOULD check row counts before and after.
* Prefer columnar or vectorized operations for bulk transformations; do not force
  vectorization when row-oriented logic is clearer or required by an external API.
* **Aggregates are derived, never accumulated.** Keep the most granular (atomic) table as the
  source of truth. Re-derive summaries from it after any filter; never obtain them by combining
  previously aggregated rows.
* **Entity keys are stable.** Each entity has a documented key (`PROJECT_CONTEXT.md` §5) that is
  identical across modules. Any rule that merges or splits entities is a named constant and is
  applied in one place.
* Use Parquet for intermediate and output tables (it preserves dtypes). Use TSV/CSV only for
  human-facing exports or tool interoperability.

## 8. Type hints

* Add type hints to public functions and non-trivial internal functions.
* Use built-in generics (`list[str]`, `dict[str, int]`) and `X | None`, not `typing.List` or
  `Optional[X]`.
* Do not add annotations that reduce readability or that cannot express the type meaningfully.

## 9. Documentation

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
* Comments explain *why*, not *what*. Use them for non-obvious decisions, assumptions,
  performance considerations and tricky logic; do not comment self-explanatory code.
* Keep documentation synchronized with the implementation, including `README.md` and
  `PROJECT_CONTEXT.md` when an interface or design decision changes.

## 10. Robustness

* Validate important inputs and fail early with informative error messages that name the
  offending value or column.
* Handle expected edge cases explicitly (empty inputs, missing groups, single-sample designs).
* MUST NOT silently swallow errors. Avoid broad `except Exception` unless there is a specific,
  documented reason, and re-raise or log with context.
* MUST NOT mutate input arguments (DataFrames, arrays, lists, dicts) in place unless this is
  intentional and documented.

## 11. Numerical and scientific code

* Make units, coordinate systems, normalization conventions and statistical assumptions explicit
  when they are not obvious.
* Avoid exact floating-point equality unless mathematically guaranteed; use tolerances.
* Handle NaN, null, infinity and empty inputs deliberately. In polars, `null` and `NaN` are
  distinct; decide which one each column uses.
* MUST NOT silently change numerical precision or relevant dtypes.
* Prefer established numerical, statistical and domain implementations over custom ones.
* **One internal coordinate convention.** Use a single internal convention (default: 0-based,
  half-open). Convert every external format at the parser, at load time, never later. Record the
  convention of each input format in `PROJECT_CONTEXT.md` §6.
* Genome build and annotation release are part of every coordinate-bearing table and are
  checked when tables are combined.
* **Every reported number states its universe**: which input set, which method or branch, and
  which filters it was computed on. A count without its universe is not a result.
* **Classifications carry their rule.** A column produced by a documented classification rule
  (as opposed to a direct observation) MUST be reported together with the identifier or version
  of the rule that produced it.
* Thresholds and analysis parameters are named constants in one configuration module, each with
  a short comment giving its rationale or source. MUST NOT hard-code them inline.
* Statistical practice (test choice, multiple testing, effect sizes, reporting) follows
  `STATISTICS_GUIDELINES.md` when present.

## 12. External tools and subprocesses

* Call external tools with `subprocess.run([...], check=True)` using an argument list.
  MUST NOT use `shell=True` with interpolated values.
* Check that required executables exist (`shutil.which`) before starting long work, and fail
  with a message naming the missing tool.
* Log the exact command line and the tool version for every external call that produces a
  result; store the version with the outputs (§13).
* Capture stderr and include its tail in the error message when a tool fails.
* Parse tool outputs in dedicated parser functions that convert to internal conventions (§11)
  and return typed tables.

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

## 14. Configuration and paths

* Library code MUST NOT contain absolute or machine-specific paths. Paths enter through
  function parameters, command-line arguments or configuration.
* Resolve the project root from an environment variable (e.g. `PROJECT_ROOT`) or by locating
  `pyproject.toml` upwards; keep machine-specific settings in an untracked local config file.
* Command-line drivers use `argparse` (or `typer` if already a dependency), provide `--help`, and
  expose every analysis parameter that a user may need to change.
* MUST NOT store secrets, tokens or credentials in code, configuration files or notebooks.

## 15. Reproducibility

* Expose random seeds as parameters and use local generators
  (`numpy.random.default_rng(seed)`); MUST NOT rely on implicit global random state.
* Results MUST NOT depend on the number of workers, parallel execution order, filesystem
  listing order or iteration order of unordered collections. Sort explicitly where order
  matters; derive per-task seeds from the master seed and a stable task key.
* Document operations whose results may vary across runs or environments.
* Keep the environment lockfile (§3) current.

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
* Keep tests independent: no shared mutable state and no reliance on execution order.
* For pipelines, SHOULD include an end-to-end test on a minimal fixture whose expected output is
  checked in.

## 17. Style

* Follow PEP 8. Formatting is handled by `ruff format`; do not hand-format against it.
* Maximum line length: 99, configured in `pyproject.toml`.
* Keep imports at module level unless a lazy or conditional import has a concrete reason.
  Group them as standard library, third-party, local (ruff rule `I`).
* Use f-strings for interpolation, except in `logging` calls, which use lazy `%` arguments.
* SHOULD configure `pre-commit` with ruff (lint and format) and `nbstripout`.

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
* Prefer pure functions; isolate I/O and side effects at the edges (scripts, loaders, writers).

## 21. Compatibility and interfaces

* Treat public APIs, documented output schemas, filenames, column names and serialized formats
  as interfaces. List them in `PROJECT_CONTEXT.md` §7.
* MUST NOT rename, remove or change the semantics of an interface unless the task requires it.
  When it must change, state the compatibility implications explicitly and record the change in
  the decisions log.
* A change to a rule that defines an output (filter, threshold, classification rule, key) makes
  previous outputs stale: say so and record it in the decisions log.
* Code shared across repositories is imported from its home package, never copied.

## 22. Notebooks

* Notebooks orchestrate and plot. Logic lives in the project package and is imported.
* The first code cell sets up paths (from `PROJECT_ROOT` or the local config, §14) and enables
  `%load_ext autoreload` / `%autoreload 2`.
* Figures go to `results/figures/<notebook-prefix>_<name>.svg` (vector by default) and follow
  `FIGURE_GUIDELINES.md` when present.
* Notebooks MUST run top to bottom on a fresh kernel. Commit them with outputs stripped.
* Notebooks are excluded from ruff and exempt from the docstring rules.

## 23. Version control

* Keep commits focused on one logical change, with an imperative summary line
  (`Fix off-by-one in GFF parser`).
* Maintain a `.gitignore` covering data, results, environment directories, caches and local
  config.
* MUST NOT rewrite published history or force-push shared branches unless the user asks.
* Message commit MUST be added at the end of a meaningful operation or code modification.
* The user MUST accept the commit first, based on the diffs and commit message suggested.

## 24. Validation

* After modifying code, run the most relevant tests and static checks (`pytest`, `ruff check`,
  `ruff format --check`): targeted tests first, the broader suite when appropriate.
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
* Do not create documents, summaries or example scripts beyond the documentation set of §2
  unless asked.
* When presenting edits to an existing file, show only the changed sections unless the file is
  short.

## 26. Related guidelines

Available in the same repository; copy the ones relevant to the project:

* `STATISTICS_GUIDELINES.md` — experimental design, test choice, multiple testing, reporting.
* `FIGURE_GUIDELINES.md` — chart choice, colour, layout, export.
* `RNASEQ_GUIDELINES.md` — RNA-seq processing and differential expression.
