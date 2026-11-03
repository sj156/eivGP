# Application data collections

Keep each related collection together:

- `ocean_silicate/n200/`, `n500/`, `n1000/`: cohort files, folds, masks and manifests.
- `adni/` and `oai/`: data contracts and private prepared fold collections.

Files within each cohort folder use their descriptive names directly. Scripts
and helpers remain together in `tools/`; they do not need one folder per script.
External prepared data directories remain supported.
