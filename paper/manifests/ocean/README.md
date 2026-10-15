# Ocean inputs, workflow identities and validation

- `upstream.json`: source provenance and current workflow defaults.
- `inputs.csv`: counts and hashes for frozen 200/500 inputs and the raw 1000 cohort.
- `validation.json`: current checks plus historical validation records.

Version 0.4.0.9008 fixes the paper protocol at four chains, 50,000 initial
iterations per chain including 5,000 burn-in. Pilot is 4/80/30; quick is 1/80/30.
Initial sampling knobs are not CLI options. Cores remains a shared fitting budget;
analysis and HTML use one worker. Explicit continuation preserves the initial
chains and may select an added length through --continuation-iter.

Default outputs start at ~/eivGP-ocean/<dataset>/<profile>/. Compatible initial
fits reuse checkpoints. Incompatible or untracked old results are preserved and
initial fitting creates a new _run_<timestamp> directory, printing its path and
an explicit --output argument. A marker next to the default root records the
latest selected run; default continue/analyze/report follow it. Explicit --output
selects an exact path. Analysis and continuation never create replacement runs.
The low-level identity/receipt safeguards remain in place.

Current validation: 954 package assertions pass; both literal main-README
profile-only short-chain commands pass (two folds/O20, eight initial tasks,
two +40 continuations, 8-core peak, serial analysis, automatic HTML). A separate
real short-chain test preserves old files, creates a new result root and follows
that root for default continuation (+5; pooled draws200 to220). Paper500 settings
are plan-checked, and sampling CLI overrides are rejected. No formal50k fit,
RAM scheduling or Linux runtime test is included. The complete versioned HTML
is rebuilt from existing short-chain results with fitting calls guarded.

Sampler, non-Ocean function definitions, synthetic studies and data are unchanged.
The raw1000 cohort is prepared separately using the R fold script. This revision integrates the verified private workflow into the public repository.
Historical validation records retain the scope and status of earlier checks.
