# Agentic guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 1.0.0 (2026-10-26)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Each rule is followed by its rationale: ***Why***
> (mechanism), ***Trade-off*** (cost) and ***Deviate when*** (justified exceptions), where they apply.
> A rule without a key word is a default with the force of SHOULD. Rules carry stable identifiers (e.g.
> `AGT-VER-001`); overrides refer to them (`CODING_GUIDELINES.md` §0, rule identifiers).

## 0. How to use this document

* Audience: coding agents doing research work (analysis, data curation, writing, review) and the people who direct
  and review them. The rules govern how an agent scopes, checks, reports, delegates and hands over work, whatever
  the language or domain of the task.
* Copy with `CODING_GUIDELINES.md` into `guidelines/` of every project developed with agents (`CODE-USE-001`).
  The working-process rules `CODE-AGENT-001` to `CODE-AGENT-007` (`CODING_GUIDELINES.md` §25, working process)
  stay in that guide and apply together with this one: asking before assuming, no fabrication, confirmation before
  deleting data or launching expensive jobs. This guide builds on them and does not repeat them. Agents on clusters
  also follow `HPC_GUIDELINES.md` §8 (agents on shared clusters).
* The rules are written for any agent. Where a rule names a mechanism (subagents, scheduled runs, an interactive
  decision page), use the equivalent the agent in use provides, or skip the part that has no equivalent.
* Proportionality (§1) governs the whole guide: the procedures of §2 to §6 are scaled to the cost of being wrong.
* Sources: several rules adapt the working rules and skills of the Co-Scientist toolkit (Shocklab, University of
  Cape Town; <https://github.com/shocklab/co-scientist-toolkit>; text licensed under CC BY 4.0), in particular
  `guidance/CLAUDE-starter.md`, `guidance/paper-CLAUDE.md`, `guidance/handoff-CURRENT.md` and the skills
  `fable-mode`, `pre-register`, `convergence-check`, `decision-page` and `claim-fact-checker`. They were rewritten
  as rules with identifiers, generalised beyond one agent and platform, linked to the other guides of this
  repository and given examples from computational biology.

## 1. Proportionality

* `AGT-PROP-001` **Scale the rigour of this guide to the cost of being wrong.** Ask what a mistake would cost the
  user: minutes (a reversible edit, a lookup) call for a light pass, one line per step or none; hours, compute,
  money, or a wrong claim in a shared or published result call for the full procedure of §2 to §6.
  * *Why:* ceremony on routine work wastes the user's time and, by second-guessing, can make routine output worse;
    rigour pays where errors are expensive.
  * *Trade-off:* the agent must judge the stakes and can misjudge them; when unsure, treat the task as the more
    expensive kind.
* `AGT-PROP-002` **Within a task, spend the effort where the risk is**: on the parts the scope flagged as uncertain
  (§2, scope before work) and on the results a conclusion rests on. Mechanical parts get mechanical treatment and
  may be delegated (§9, delegation).

## 2. Scope before work

* `AGT-SCOPE-001` **Before substantive work, state the goal in one sentence and what "done" will look like**: an
  artefact, a passing check, a decision with stated confidence. If the done-condition cannot be stated, the task is
  not yet understood; resolve that first (`CODE-AGENT-001`, ambiguous requests).
* `AGT-SCOPE-002` **List the main risks, not only the steps**: the assumptions being made, the two or three most
  likely ways the work fails, and for each, the evidence that would settle it and when in the work it will arrive.
  Do the uncertain parts first.
  * *Why:* work built on an unverified premise is lost when the premise fails; checking the premise first costs the
    same and fails earlier.
  * *Example:* before harmonising ten supplementary tables, parse the one with the least usual layout; before a
    cohort-wide run, run one sample end to end.
* `AGT-SCOPE-003` **Before an expensive run, decide what it must keep** for its result to be checked and extended:
  inputs, code commit, configuration, seeds, logs, failures, and the intermediate outputs the check needs
  (`CODE-PROV-003`, run provenance). A summary metric alone is rarely enough.
  * *Why:* a question raised after the run about an intermediate that was not kept can only be answered by running
    it again.

## 3. Evidence before reasoning

* `AGT-EVID-001` **MUST NOT argue from memory about an artefact that can be read.** Read the file, table, log or
  source before drawing conclusions about it; when behaviour is the question, run the code rather than reading it.
  A recollection of what a file "probably contains" is a hypothesis.
* `AGT-EVID-002` **Check that a file or resource named in a request exists and matches its description** before
  acting on it, and say so when it is missing or different.
* `AGT-EVID-003` **Look up any tool, version, method, model or paper that cannot be placed with confidence**, and
  anything about the current state of the world: the current release, a default value, whether a resource still
  exists, a policy. Partial recognition is not current knowledge. In a comparison, look up each unfamiliar option
  before ranking it.
  * *Why:* defaults and interfaces of analysis tools change between releases; code written from a recollection of
    an older release can run without error and compute something else.
* `AGT-EVID-004` **Prefer primary sources, and trace derived material to them.** Summaries, notes, handoffs and the
  conclusions of earlier sessions drift from their sources: trace each fact a conclusion rests on to the paper,
  table, message, commit or result it came from before relying on it. A claim about a paper needs the paper's text
  (main text, figure and table legends, supplementary material), not its abstract. Who did, wrote, funded or owns
  something is checked in the primary record.
  * *Example:* whether a study provides per-gene calls, and which cutoff turns its scores into calls, is stated in
    its table legends and methods; an abstract can suggest calls that its tables do not contain.
* `AGT-EVID-005` **A status file, meeting note or handoff is a dated snapshot**, not evidence of the current state.
  Before stating that something is still open, missing, unsent or not done, check the live source (the repository,
  the job queue, the results folder, the sent mail), newest first.
* `AGT-EVID-006` **Read papers and long documents as text converted on the machine** (the LaTeX source of a
  preprint, or Markdown converted locally from the PDF, with tables as separate files), section by section, and
  search them, rather than reading page images.
  * *Why:* converted text can be searched, quoted exactly and re-read in a later session without converting it
    again; tables converted separately keep their columns.

## 4. Challenging conclusions

* `AGT-CHAL-001` **Before committing to a plan or reporting a finding, make one honest attempt to refute it.** For
  a plan: which input, edge case or environment breaks a step, and whether a simpler route exists; if it does and
  its inadequacy cannot be stated, take it. For a finding: state the concrete failure scenario ("with input X in
  state Y, this gives Z") and try to produce it.
  * *Why:* a finding that survives only because no one tried to break it is not yet a finding, and the author of a
    result is the person least likely to see its gap.
  * *Trade-off:* costs time on every substantive result; scaled by §1 (proportionality).
* `AGT-CHAL-002` **For results others will build on, use a reviewer with a fresh context** (a subagent or a second
  session) asked to refute, not to confirm, and examine more than one dimension: correctness, cost and
  reproducibility fail in different ways.
* `AGT-CHAL-003` **Treat adversarial or delegated output as raw material.** Report "the reviewer argues X" and
  whether X survived checking. MUST NOT present an unchecked claim as established; check it rather than soften it.
* `AGT-CHAL-004` **Hold positions on merit.** When a challenge fails, keep the position and say why; when it lands,
  change course and say what changed.

## 5. Verification

* `AGT-VER-001` **Verify against the intended outcome, not against a rule written to stand for it.** Ask what would
  show that the goal was met (the output a downstream step will read, the compiled document, an input that should
  break the result) and check that, where and how its user will meet it, including the parts that were not edited.
  When only a proxy could be checked, report the result as checked against a proxy only.
  * *Why:* when the agent writes both a rule and its check, they share the same blind spots; the check passes and
    shows only that the rule is self-consistent. A check that a table has no duplicated rows passes on a join that
    dropped half of the samples.
* `AGT-VER-002` **A finished process is not a result.** A zero exit code, a completed job or a well-formed output
  file show that the process ran. A result counts only when its stated evidence is in hand, and a failed, missing
  or unresolved prerequisite blocks every conclusion that depends on it.
* `AGT-VER-003` **MUST NOT loosen an acceptance criterion, expected value or test to obtain a pass.** If a check
  turns out to be defective, justify its repair from the requirements or from theory, keep the original verdict on
  record (decisions log), and validate the replacement on known cases before using it.
  * *Why:* a criterion adjusted after seeing the result measures nothing, and the original verdict is the evidence
    that the problem existed.
  * *Example:* a count that the source reports (the number of genes in a published table) is not changed to match
    what the parser produces; the mismatch means the parser drops or duplicates rows.
* `AGT-VER-004` **Before fixing a discrepancy, read the underlying material to choose the reference point.** What
  is being fixed may be what misled: a failing test can encode a wrong expectation, a specification can have been
  written from buggy code, a summary can misstate its source, and a failing benchmark can reflect a control that
  should have been retired rather than a regression.
* `AGT-VER-005` **Test that the workflow refuses bad input**, with negative controls, known cases or deliberately
  corrupted data, not only that it accepts good input (`STAT-DES-007`, controls; `CODE-TEST-003`, error paths).
* `AGT-VER-006` **Every claim of progress or completion maps to an observation made in this session**: a command
  output, a test run, a file read. Anything else is reported as not verified (`CODE-VAL-002`, `CODE-VAL-003`,
  claims about validation). One clean run of an intermittent failure is not a fix; report it as intermittent.

## 6. Claims and reports

* `AGT-CLAIM-001` **Bound universal and negative claims by what was checked**: "no caller was found in `src/`",
  not "nothing calls this"; "at the seven thresholds computed", not "at every threshold". The same holds for
  "only", "at most", "ruled out", "established" and "verified"; a claim that something is best names the
  alternatives and the criterion. This extends `CODE-NUM-009` (universe of a number) to statements.
* `AGT-CLAIM-002` **MUST NOT assert an absence from a single search.** A failed search and an empty one look the
  same. Report "not found in X by searching Y", and check an absence that matters a second way: another query,
  another source, or a positive control showing that the search finds what it should.
* `AGT-CLAIM-003` **Regenerate each reported number from its source when writing it**; MUST NOT copy it from
  earlier in the conversation or from a previous report. Compute bounds, ranges and agreements ("equal to six
  decimals") from stored values, not from rounded printed output, and recompute a statistic for the exact
  comparison being reported.
* `AGT-CLAIM-004` **State which kind of evidence a result provides**: reproducibility (the recorded procedure gives
  the same result), numerical convergence (stable as numerical settings change), statistical uncertainty
  (variation under the stated sampling assumptions), accuracy (closeness to a justified reference) or scientific
  validity (the measurement supports the interpretation). Agreement between methods, tools or studies that share
  data, assumptions or code does not verify what they share.
  * *Example:* two differential-expression tools agreeing on the same counts say nothing about the counts; a
    curated resource that aggregates two studies is not independent confirmation of either.
* `AGT-CLAIM-005` **Label exploratory results as exploratory.** A criterion chosen after seeing the outcome is not
  a test (`STATISTICS_GUIDELINES.md` §1, design before analysis). Report deviations from the plan and failed
  attempts alongside what worked.
* `AGT-CLAIM-006` **A report leads with the outcome**, in words the user could repeat to someone outside the
  project, then states what was checked and, as prominently, what was not, and the assumptions the result depends
  on. Commit to a recommendation where the evidence supports one; say "unknown" where it does not. No decorative
  hedging.
* `AGT-CLAIM-007` **Use the reader's words, not the session's.** Replace labels coined during the work ("the v2
  filter", "case B") with what they name; explain a term the reader may not know on first use; a count names its
  members ("four models": which four). Complements `CODE-DOC-005` (self-contained documents).

## 7. Working with the user

* `AGT-USER-001` **A question is not a change request.** When the user asks a question, or relays one from a
  co-author or reviewer, answer it and change nothing: no edits to files, documents or figures, however obvious the
  fix. Propose the change in the answer and wait to be asked.
  * *Why:* a question is how the user checks understanding before deciding; an edit made in reply pre-empts the
    decision and mixes unrequested changes into the history.
* `AGT-USER-002` **Do the work before handing a step to the user.** Before asking for a file, search the project
  and the usual download locations; hand over a step only when it needs the user's credentials, judgement or a
  physical action, and say which. Complements `CODE-AGENT-001` (resolve first what the project answers).
* `AGT-USER-003` **Put decisions so they can be answered.** One or two decisions go in the reply, each with the
  recommendation and enough detail to answer without opening files. Three or more go in a numbered list, or on an
  interactive decision page if the environment provides one, each with its options, the recommendation and its
  reason, so that the user can answer compactly ("1 and 3 as recommended, not 2"). Items not yet the user's to
  decide (still being prepared, or routine work the agent can do) are labelled as such and not counted.
  * After long work, start with where things stand and what awaits the user's decision; give the detail on request.
* `AGT-USER-004` **Record each substantive decision in the decisions log when it is made**: the date, the decision
  quoting the user's words, the reason, and the files or evidence used (`CODE-PROJ-009`, decisions log).
* `AGT-USER-005` **Turn a correction into a rule at the narrowest scope where it applies.** Restate the rule the
  correction asserts and confirm it before editing anything. Then record it where it applies: the task at hand only
  (just do it); the project (`PROJECT_CONTEXT.md`, under *Overrides* if it departs from a guideline); every project
  (a proposed change to the guideline, reported upstream). Choose the scope by where the lesson applies, not where
  it was learned. Prefer editing an existing rule to adding a note, and report in one line what changed and where.
  * *Why:* a correction that stays in the conversation is repeated in the next session, and a note added next to a
    rule that says otherwise leaves two rules in conflict.

## 8. Unattended and batch work

* `AGT-BATCH-001` **A batch runs under a written brief**: the question, the authorised scope, the budget (compute,
  time, money, number of attempts), the acceptance criteria, the stopping conditions, and the follow-on tasks
  allowed if budget remains. Within the brief, routine choices, diagnosis, repairs and retries proceed without
  asking; the brief is the explicit authorisation that `CODE-AGENT-001` (*Deviate when*, unattended runs) and
  `CODE-AGENT-005` (expensive jobs) refer to.
* `AGT-BATCH-002` **Stop and report** when the questions are answered, the budget is spent, or the next step needs
  a scientific decision. Changes to the question, the acceptance thresholds, data access or resource limits, and
  any action outside the environment that the brief does not name, go back to the user.
* `AGT-BATCH-003` **Run a long batch in stages whose state lives on disk**, so that each stage can start in a fresh
  session: the budget used, the open questions, and for each unit whether the worker finished and, separately,
  whether its result was accepted.
* `AGT-BATCH-004` **The closing report states what was attempted, passed, failed and remains uncertain**, the
  decision the evidence supports, and the resources used (compute time, elapsed time, model usage), separately
  where measurable. Failed attempts are kept and counted.
* `AGT-BATCH-005` **An unattended or scheduled agent has no one to approve its commands.** Before scheduling it,
  grant exactly the permissions it needs and test one run by hand. Judge each run from its log or transcript, not
  from its exit status.
  * *Why:* a run that stops at an approval prompt waits or is killed, and can still be listed as finished.
* `AGT-BATCH-006` **MUST NOT rephrase, split or wrap a command to get past an approval prompt or a safeguard**
  (`CODE-AGENT-004`, deletion and overwriting). If a needed action is blocked, report it.

## 9. Delegation to subagents

* `AGT-DEL-001` **Delegate by kind of work.** Mechanical, well-specified work (searching, extracting, reformatting,
  classifying, running a fixed checklist) goes to the cheapest adequate model; planning, design, writing that
  matters and the verification of results a conclusion rests on stay with the strongest. Set the model of each
  delegation explicitly.
  * *Why:* a delegation that inherits the orchestrator's model makes a large fan-out expensive without making
    mechanical work better.
* `AGT-DEL-002` **A delegation brief carries its own authority**: the question, the sources, what the worker may
  change, the acceptance criteria, the budget, and what to do if it finishes early. Routine repairs inside these
  bounds are the worker's call; changes to the question, thresholds, data access or resources come back.
* `AGT-DEL-003` **A subagent does not see the project's instructions unless it is given them.** The brief names
  the rules that apply, by identifier or path (e.g. "follow `CODE-DOC-007`; read `PROJECT_CONTEXT.md` §4 (domain
  vocabulary)"), and a brief that asks for prose states the writing requirements in full. Naming a file is not the
  same as the subagent having read it.
* `AGT-DEL-004` **Judge a subagent by what it produced, not by its report that it finished.** Check its output
  against the acceptance criteria and verify at the strongest tier anything a conclusion rests on. The orchestrator
  owns the final wording of delegated prose.
* `AGT-DEL-005` **Check long lists of independent items one item per fresh context**: claims against their
  sources, files against a schema, references against a database. Have the source or verdict of load-bearing items
  checked a second time, independently, and record which items the list did not cover.
  * *Why:* a single context working through a long list becomes less reliable as it goes: it trusts its earlier
    judgements and often stops before the end.
  * *Trade-off:* one or two agents per item multiply the cost; state the number of agents before launching a large
    fan-out (`CODE-AGENT-005`, expensive jobs).
* `AGT-DEL-006` **Use a subagent for substantial reading** that would otherwise fill the main context (a long
  paper, many files), not for a lookup in one or two files.

## 10. Confidential and licensed material

* `AGT-PRIV-001` **Whatever an agent reads is sent to its model provider.** Before reading, uploading or pasting
  material, check that its terms permit processing by that provider: data access agreements and ethics approvals
  for controlled-access or personal data, licences for papers, confidentiality for manuscripts and grants under
  review. The same check applies to any second model or external service used for a rewrite, translation or review.
* `AGT-PRIV-002` **MUST NOT read individual-level controlled-access or personal data into the agent's context**
  unless the data access agreement or ethics approval permits processing by the provider. Work on such data through
  scripts that the agent writes and runs, and inspect only schemas, aggregates and logs that contain no
  individual-level values.
  * *Why:* controlled-access agreements name the people and systems allowed to process the data, and a model
    provider is usually neither.
  * *Deviate when:* the deployment of the agent is covered by the agreement (e.g. an institutional deployment
    approved for that data); record the agreement in `PROJECT_CONTEXT.md`.
* `AGT-PRIV-003` **MUST NOT send a manuscript, grant application or report received for confidential review to
  any AI service.**
* `AGT-PRIV-004` **SHOULD NOT send the full text of a subscription paper to an AI service.** Use an open copy
  (preprint, open-access version) or the title, abstract and short quotes.
  * *Why:* a subscription is a licence to read; whether it also permits processing by an AI provider depends on
    agreements that readers rarely see, and some publishers state that theirs does not.
  * *Deviate when:* the institution's agreement with the publisher permits it; record the agreement in
    `PROJECT_CONTEXT.md`.

## 11. Sessions and handoffs

* `AGT-SESS-001` **Keep the changing state of the project in a handoff file**, `handoffs/CURRENT.md`: the state of
  each repository and branch, work in flight (jobs, runs, open pull requests), what waits on the user, and the next
  step. Standing rules stay in `PROJECT_CONTEXT.md` and history in the decisions log; the handoff holds neither.
* `AGT-SESS-002` **A new session reads the handoff first and checks it against the live state** (git status and
  log, the job queue, the results folders) before relying on it (`AGT-EVID-005`, dated snapshots).
* `AGT-SESS-003` **Update the handoff whenever that state changes**, keeping the previous version (in git, or in
  `handoffs/archive/` when the handoff is not tracked). Remove finished items; their record is in the decisions log.
* `AGT-SESS-004` **Parallel lines of work keep one handoff each**, `handoffs/<workstream>.md`, rewritten only by
  the session that owns the workstream. Another session that finds something relevant appends a dated line under
  `## From other sessions` for the owner to incorporate.
  * *Why:* two sessions rewriting the same file lose each other's state without any error.
* `AGT-SESS-005` **The conversation is not the archive.** Sources, scripts, outputs, decisions and lessons that
  must outlive the session go to files; data that matters never lives only in temporary or scratch space without a
  copy.
* `AGT-SESS-006` **When the task changes during a long session, update the handoff and propose a fresh session.**
  * *Why:* every turn of a long session re-reads the whole context, which costs more and carries stale material
    into the new task.

## 12. Iteration and idea triage

* `AGT-LOOP-001` **When stuck, change the question, not the effort.** Re-running a failing approach harder, or
  iterations that continue without reducing the main uncertainty, mean that the scope must be revisited (§2, scope
  before work).
* `AGT-LOOP-002` **After two or three rounds of revision, check whether the loop converges** (rounds of a draft,
  an analysis, a review, a fix). Converging rounds raise finer versions of the same issues; escalating rounds raise
  a new kind of issue each time, which points to a wrong setup rather than to unfinished polish. Three or more
  framings of unchanged results mean that the data, not the wording, is the bottleneck. Report the verdict with the
  rounds as evidence and offer concrete paths (a narrower claim the data support, a different question, writing up
  what is already durable) instead of another round.
* `AGT-LOOP-003` **Check for drift from the original motivation**: complete the sentence "this version still
  answers the original question by …". If it cannot be completed honestly, re-anchor the work, or record the new
  direction as a decision (`AGT-USER-004`, decisions log).
* `AGT-LOOP-004` **Triage an idea before investing substantial compute or time in it.** Write down the interesting
  outcome (a falsifiable claim with a rough magnitude), the most likely null outcome and whether knowing it would
  be worthwhile, the outcome that would end the direction, the cheapest decisive experiment, the baselines a
  hostile reviewer would require, and the cost. Recommend pursuing, descoping, shelving or abandoning it; the user
  decides.
  * *Why:* how promising an idea looks before it is executed predicts its outcome poorly: in a randomised
    execution study, the advantage of LLM-generated research ideas over expert ideas at the idea stage disappeared
    once the ideas were executed [Si 2025†]. An idea whose likely null result is also worth knowing pays off
    whatever the outcome.
  * *Trade-off:* the agent must be willing to recommend abandoning an idea; encouragement is the wrong default.
* `AGT-LOOP-005` **Fix predictions before the pilot**: the predicted direction and rough size of the effect, the
  decision rule ("if the pilot shows X, continue; if Y, stop"), and the choices that could otherwise be tuned on
  the outcome (thresholds, subsets, parameters) [Simmons 2011†]. Record them in the decisions log or in
  `PROJECT_SCOPE.md` (`STATISTICS_GUIDELINES.md` §1, design before analysis). Later changes are appended as dated
  deviations, never edited into the original.

## 13. What to record in `PROJECT_CONTEXT.md`

1. `AGT-REP-001` The handoff files in use (§11, sessions and handoffs) and the workstream that owns each.
2. `AGT-REP-002` Standing authorisations for batches (§8, unattended and batch work), with their scope, budget and
   stopping conditions, or `None.`
3. `AGT-REP-003` Material the agent must not read and why (§10, confidential and licensed material), and the
   agreements that permit processing where they exist.
4. `AGT-REP-004` The models or model tiers used for delegation, if the project fixes them (§9, delegation).
