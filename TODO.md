# To do

Work planned for the guides, recipes and tooling that is not yet in them. Each item says what is missing and where it
would go; remove an item in the change that completes it.

## Guides

- `SCRNASEQ_GUIDELINES.md`: pooled CRISPR screens read out by the transcriptome (guide assignment, perturbation
  effects, multiplicity of infection). The 4.1.0 draft mentioned them in the scope (§0) and pointed to §15, but no
  section exists; the mention was removed in 4.1.0 until the section is written.
- Established pipelines, next round (`CODING_GUIDELINES.md` §12, *Established pipelines*): whether the levels of
  `CODE-PIPE-002` extend to analysis libraries (statistics, enrichment) and when exploratory code is promoted to a
  pipeline step; a pipeline baseline for analyses after the processing stages and for figures; criteria for an
  unreleased pipeline beyond a pinned commit with passing tests.
- `SCRNASEQ_GUIDELINES.md` `SC-QUANT-001`: motivating case for the pinned fork of `nf-core/scrnaseq` (nuclei with
  `GeneFull_Ex50pAS`), to be written from the first project that runs it.

## Tooling

- `tools/check_deviations.py`: read the `pipeline_deviations` block of `PROJECT_CONTEXT.md` (`CODE-PIPE-004`), check
  each typed return condition (`schema_accepts` against the `nextflow_schema.json` of a release, `first_release`,
  `pipeline_includes`, `issue_closed`) and list the deviations whose condition is met (`CODE-PIPE-005`); `manual`
  conditions listed for review.
