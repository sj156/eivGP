# OAI reproduction and storage guide

The application uses the same flat tools layout, prepared-fold contract, package
samplers, numerical prediction functions and resource allocation as ADNI/Ocean.
OAI-specific labels and display adjustments live in `tools/oai-tools.R`. That
adapter reuses the existing prepared-fold engine in `tools/adni-tools.R`; it does
not duplicate fitting or scoring implementations. Reusable numerical and parallel
functions remain in the package.

## 1. Install and supply inputs

Follow the repository installation instructions, reinstalling this checkout.
Place the privately supplied 46 files directly in `data/oai/`; see
[the data contract](../../data/oai/README.md). An external directory works with
`--data-root=/local/oai`. Paths in the code are relative to the checkout; commands
can be invoked from another working directory.

The default scientific design uses baseline quantitative MRI cartilage thickness
as U, fixed-location medial radiographic joint-space width as Y, continuous
age/BMI/height as X, and eight native central-reading radiographic grades as C.
Only the documented adjacent lateral-femoral-sclerosis grades are pooled in the
supplied cohort. Metadata provides all coordinate labels, grade maps and
selection provenance; another supplied p/q/d configuration does not require
rewriting the fitting code.

## 2. Paper execution

```sh
Rscript tools/oai-eivGP.R --plan=true --profile=paper --cores=16
Rscript tools/oai-eivGP.R --profile=paper --cores=16
Rscript tools/oai-competitors.R --profile=paper --cores=16
Rscript tools/oai-report.R --profile=paper --cores=16
```

Run the fitting groups sequentially against the same output root. `--cores` is
the sole resource control and is capped at detected capacity. New fits use
min(cores,4) chains; above four cores, multiples of four permit complete resident
EIV groups. Each native competitor job uses one core. Lower the budget if memory
is limited. Changing only resources for saved-only reporting does not refit.

Paper settings initially use 50,000 total iterations per chain, burn 5,000,
HMC O0/O20 and ESS O50. Prediction selects 250 retained states per chain.
`--ess-min=100` optionally continues the saved states until monitored scientific
bulk and tail ESS exceed 100; it does not certify convergence. Explicit `--iter`,
`--burn`, `--stages` and sampler overrides are frozen with the run identity.

## 3. Pilot and quick checks

`--profile=pilot` uses the same fivefold workflow with a single 1,500-iteration
stage, burn 500 and 250 prediction states per chain. With at least four cores this
is four chains retaining 1,000 states each. Competitor controls remain the regular
paper controls. Paper mode is a separate output root, not a reinterpretation of
historical pilot fits.

`--profile=quick` shortens the synthetic check and competitor optimization; it is
not a scientific performance comparison. A local synthetic check exercises the
actual EIV, competitor and report CLIs, including saved-fit reuse:

```sh
Rscript tests/application/check-oai.R --cores=4
```

No OAI records are used by that check.

## 4. Saved-only rendering and outputs

```sh
Rscript tools/oai-report.R --profile=paper --cores=4
Rscript tools/oai-eivGP.R --profile=paper --cores=4 --reuse-only=true
Rscript tools/oai-competitors.R --profile=paper --cores=4 --reuse-only=true
```

The combined saved report is `results/oai/oai.html`. Method pages are
`results/oai/reports/eivGP/index.html` and
`results/oai/reports/competitors/index.html`. The editable source is
`paper/applications/oai.Rmd`; `paper/applications/oai.html` is a data-free plan.

| Result block | Contents |
| --- | --- |
| `eiv_foldNN_O*/` | Final fit, predictions and diagnostics; earlier stages stay inside the task |
| `baseline_foldNN/`, `cc_foldNN_O*_METHOD/` | Native OI and complete-case fits and scores |
| `data/`, `_meta/`, `logs/` | Verified private input snapshots, identities, receipts and task logs |
| `reports/` | Method pages, prediction summaries, paper exports and scientific diagnostics |

The marker remains `adni-flat-v1` because OAI uses the same prepared-fold storage
engine; task names and folder shapes match ADNI. Existing unmarked result roots
retain their old paths. Use `--output` consistently for a nondefault root.

## 5. Interpretation and provenance

OI predicts all held-out participants from (X,C). CC baselines train only on
revealed calibration U and use (X,U,C); EIV conditions on the same held-out U in
f(X,U). Every CC method uses the identical supplied known-U test subset. The
workflow does not search or optimize test masks. Metadata records test revelation
separately from training calibration; outcome-selected supplied masks are
explicitly labelled exploratory in the report.

Training-U recovery uses hidden training rows and may condition on training Y;
prospective recovery uses test C without test Y. O0 has no physical-U error score.
The optional metadata field `imputation_display_fold` selects the main training-U
figure, falling back to the configured representative fold when unavailable.
Response-function and trace displays follow `--representative-fold`.

All finite methods remain in metric rankings. Asterisks identify CORE failures
or exposed flagged optimizer returns; missing diagnostics are not declared
converged. Original results and dataset-specific provenance belong in local
inputs/results, not public source files.
