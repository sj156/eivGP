# Two workflows, two notebooks

Start with a rendered example, then open the matching notebook and run it.
Every displayed result is computed by its notebook; HTML contains no manually
copied result tables. Compatible caches can be reused, with input/settings
checks. Quick runs are illustrations, not the paper's inferential results.

| Workflow | Notebook | Rendered quick example |
| --- | --- | --- |
| Both numerical studies | [01_numerical_experiments.Rmd](01_numerical_experiments.Rmd) | [HTML](../docs/replication/01_numerical_experiments.html) |
| GLODAP silicate/oxygen application | [02_ocean.Rmd](02_ocean.Rmd) | [HTML](../docs/replication/02_ocean.html) |

Both default to `profile: "quick"` and `action: "run"`. Knit the notebook
in RStudio, or from the repository root:

```r
rmarkdown::render("replication/01_numerical_experiments.Rmd")
rmarkdown::render("replication/02_ocean.Rmd")
```

The numerical notebook handles both studies (`study: "both"`) and also offers
`development`. `paper` selects the full scientific configuration within the
same notebook; do not create a separate notebook for each budget. The paper
numerical profile uses 50 datasets per setting and four chains × 20,000
iterations, including 5,000 warmup. The ocean notebook lists its application
protocol and includes the recorded continuation step in the paper workflow.

Dependencies must be installed first; see the root README. The frozen
public ocean subset is included in `data/ocean_silicate/`.

The HTML shows a short settings table, simple exploratory analysis, selected
model outputs, and relevant diagnostic/method-failure information. Routine
progress, detailed diagnostics, checkpoints, and full intermediate tables stay
in local output directories. `action: "report"` reads compatible saved results
without fitting; `plan` checks settings. Change `output_dir` or `artifact_root`
when intentionally changing a cached analysis configuration.

Preprocessing, simulation generation, reporting, and operational recovery are
supporting scripts under `codes/`. Old per-stage and historical diagnostic
notebooks are preserved in the local cleanup snapshot and Git history, rather
than presented as additional reader entry points.
