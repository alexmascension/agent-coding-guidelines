# To do

Work planned for the guides, recipes and tooling that is not yet in them. Each item says what is missing and where it
would go; remove an item in the change that completes it.

## Guides

- `SCRNASEQ_GUIDELINES.md` `SC-QC-017`: a recipe in `recipes/scrnaseq/qc.py` (per-capture Gaussian mixture, ICL
  selection, outer-tail Mahalanobis flags, per-capture report and figure) with tests, and the matching row of §16.
- `SCRNASEQ_GUIDELINES.md`: pooled CRISPR screens read out by the transcriptome (guide assignment, perturbation
  effects, multiplicity of infection). The 4.1.0 draft mentioned them in the scope (§0) and pointed to §15, but no
  section exists; the mention was removed in 4.1.0 until the section is written.
- Established pipelines, next round (`CODING_GUIDELINES.md` §12, *Established pipelines*): whether the levels of
  `CODE-PIPE-002` extend to analysis libraries (statistics, enrichment) and when exploratory code is promoted to a
  pipeline step; a pipeline baseline for analyses after the processing stages and for figures; criteria for an
  unreleased pipeline beyond a pinned commit with passing tests.
- `SCRNASEQ_GUIDELINES.md` `SC-QUANT-001`: motivating case for the pinned fork of `nf-core/scrnaseq` (nuclei with
  `GeneFull_Ex50pAS`), to be written from the first project that runs it.
- `hpc_sites/DIPC_HYPERION.md`: for each observation, the command and a minimal output (or an incident reference),
  so it can be re-verified (audit item G-H06); add the `srun` note (`SLURM_EXPORT_ENV=ALL` before `srun` steps of a
  job submitted with `--export=NONE`, `HPC-JOB-001`) once tested on the cluster.
- `CODING_GUIDELINES.md` `CODE-USE-007`: the text defines MINOR as new rules or reorganization and MAJOR as removed
  or substantially changed rules; recent rounds bumped MINOR for requirements added to a rule or scopes narrowed or
  widened while overrides still apply. Write that practice into the rule, or decide otherwise.
- References: separate the kind of evidence of each note (method paper, independent benchmark, review,
  documentation, local experience) and state the scope of comparative claims (datasets, criterion); trace each
  MUST/MUST NOT, threshold and superlative to the passage or figure of its source (audit §1.15).
- `OMICS_REFERENCES.md` [Hennig 2007]: the note is empty; no abstract was accessible (CrossRef, Europe PMC).

## Tooling

- `tools/check_deviations.py`: read the `pipeline_deviations` block of `PROJECT_CONTEXT.md` (`CODE-PIPE-004`), check
  each typed return condition (`schema_accepts` against the `nextflow_schema.json` of a release, `first_release`,
  `pipeline_includes`, `issue_closed`) and list the deviations whose condition is met (`CODE-PIPE-005`); `manual`
  conditions listed for review.
- `tools/check_guides.py`: citation-key checks also for `CODING_GUIDELINES.md`, which is excluded today (audit C22).
- `recipes/common/build_gene_groups.py`: the shipped lists predate source snapshots and cannot be rebuilt offline.
  The next networked rebuild (`python -m recipes.common.build_gene_groups --overwrite`) writes the first snapshot under
  `build/gene_groups_sources/<date and time>`; archive it outside git (the HGNC response is about 19 MB, above the pre-commit
  size limit) and record where.
