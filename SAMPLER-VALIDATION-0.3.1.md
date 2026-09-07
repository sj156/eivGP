# Sampler validation: eivGP 0.3.1

## Exactness and software checks

The prior and posterior target are unchanged from 0.3.0. The coordinate
cutoff conditional is derived in `codes/COMPUTATION_CONCORDANCE.md`.
`codes/tests/test_v031_updates.R` checks truncated-Beta inverse draws against
independent quadrature, sequential Dirichlet splits, a known binary-cutoff
conditional, tight calibration constraints, empty categories, unchanged U,
the probit likelihood-factor union, and use of the joint block-GP shortcut.
Existing checks cover dense/Schur equality and fallback, prediction covariance,
all-draw retention, and reproducible uninterrupted/continued and serial/fork
execution. Repository loader and worker-error tests also pass for both engines.

## Bounded paired computation check

One deliberately difficult existing Study I case was selected from the supplied
diagnostics: eta0_balanced, replication 2, n=100, six categories, 50 calibration
observations. The saved failure record supplied the unchanged standardized X/y,
ordinal C, calibration U, and resolved priors. Each fresh fit used seed 302050,
four forked chains, 1750 transitions, 500 warmup, and all 1250 retained draws
per chain. The old cutoff step reproduced its original retained cutoff draws
exactly. Both new fits used the same posterior and data.

| Cutoff update | U block size (subjects) | Seconds | Maximum cutoff R-hat | Minimum cutoff bulk ESS | Minimum cutoff ESS/second |
| --- | ---: | ---: | ---: | ---: | ---: |
| Previous joint ESS | 8 | 29.35 | 1.4604 | 7.85 | 0.267 |
| Coordinate Gibbs | 8 | 31.33 | 1.0091 | 340.83 | 10.877 |
| Coordinate Gibbs | 4 | 40.03 | 1.0067 | 451.30 | 11.273 |

The maximum missing-U R-hat was 1.0198, 1.0073, and 1.0054 respectively;
minimum missing-U bulk ESS was 240.94, 935.63, and 1217.51. The cutoff
transition produced the large improvement in this case. Block size 4 gave
more effective draws per sweep but only a small ESS/time advantage. Runtime
is machine/load dependent, and this is one dataset/seed, not a general mixing
guarantee or numerical-study result. No scientific-target convergence claim is
made from this check, and dictionary/other diagnostic warnings remain visible.

The default U block remains 8 subjects. Use
`sampler_control = list(u_block_size = 4L)` for a smaller-block comparison.
`threshold_update = "ess"` enables the previous cutoff transition. Both are
fitting controls, not experiment-design choices. Compare settings through fresh
fits; continuation freezes the chosen controls and requires sampler 0.3.1.

No manuscript experiment settings, frozen datasets, existing results, or
diagnostic thresholds were changed. The installed package contains reusable
model APIs and tests, not this case-specific benchmark or experiment runners.
