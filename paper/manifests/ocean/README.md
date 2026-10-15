# Ocean input and workflow identities

- `upstream.json`: source repository/commit and Ocean integration scope.
- `inputs.csv`: row counts and SHA256 for 33 input CSVs: 16 each for the 200/500 frozen designs and one raw 1000-profile cohort.
- `validation.json`: checks executed for this integration and their exact scope.

The package's installed `module-manifest.csv` records source module hashes. Output roots additionally store package/functions/data identities, per-task configs, checksums, completion receipts and fitting settings. The example HTML uses a fresh, complete 200-profile five-fold pilot run, including continuation, all diagnostics and Supplementary figures. Short-chain numerical flags remain visible. The updated 500 and raw 1000 inputs have been checked; full-budget fits remain separate. The standalone fold script writes generated data under `reproduction/ocean/prepared/`, outside the versioned raw cohort. Original fitting archives are preserved outside this repository.


The development branch merged upstream `ebc297b` before aligning the
Ocean entrypoints. Experiments now run with `Rscript paper/applications/ocean.R
--action=run --profile=paper --dataset=ocean500 --cores=4`; the Rmd defaults to
`action=report` and rejects fitting. CLI defaults use
`results/ocean/<dataset>/<profile>/`. A fresh five-fold, short-chain CLI run
completed 65 fits/tasks. A CLI empty-path bug found during postprocessing was
fixed and regression-tested; the corrected full command reused all 65 tasks.
The separate report render produced 19 tables, 10 embedded images and nine
Supplementary blocks while fitting, continuation and comparator fitting were
guarded against execution. Saved fit/receipt hashes and modification times were
unchanged. All 791 package assertions passed. The three implementation modules
match the merged source; no Ocean fitting or sampler method was changed by the
entrypoint alignment. These checks were completed on the development branch before the approved
public integration.

The main commands, CLI defaults and Rmd dataset parameter now use `ocean500`.
Its no-argument plan was verified against 500 records and 65 tasks. The versioned
HTML continues to illustrate the separately verified 200-profile short-chain
run; no 500-profile formal fitting was performed in this revision.

The public checkout was separately installed and checked before publication:
791 package assertions, the default 500-profile/65-task plan, and a guarded
report-only render of the complete saved 200-profile pilot all pass. The report
contains 19 tables, 10 embedded figures and nine Supplementary blocks. Existing
fits and receipts remain unchanged.


Stage separation (0.4.0.9004): `fit` saves checkpoints and basic package
diagnostics; `analyze` reads saved chains using one worker; `continue` explicitly
requests the existing one-round, screen-selected extension. `run` aliases `fit`.
The diagnostics-worker option and Ocean diagnostic forking were removed.
The main README includes literal Terminal commands for a one-fold O20 check
using the 200-profile cohort and three CC comparators. Formal commands remain
configured for ocean500. The complete example HTML is rebuilt from existing
short-chain fits without new full-five-fold MCMC. This staged revision is published following explicit user approval.

This staged public checkout was independently installed and tested: 791 assertions
pass, and all six literal main-README Terminal commands complete successfully
on ocean200/fold1/O20 with the three CC comparators. Fit/analysis receipt hashes
verify after the explicit extension and HTML render. No formal 500-profile fit
was launched for this publication check.


Configurable continuation (0.4.0.9005): `--continuation-iter` sets added
iterations per chain when running `continue`; defaults remain paper 20,000 and
pilot 80. Each extension stores its budget in `continuation_settings.rds`, which
is included in the hashed task receipt. Initial-fit identity is independent of
this later choice. Analysis and reports read saved lengths, including legacy
fixed-budget results. A conflicting request for an existing extension is rejected.
The six literal main-README commands pass with a 40-iteration extension: four
chains add 160 pooled draws (200 to 360). Default reuse, conflict rejection,
serial reanalysis and HTML are verified. The complete saved-result example
retains 19 tables, 10 images and nine Supplementary blocks. No new full-budget
or five-fold MCMC run was performed; this update was explicitly approved for public publication.
