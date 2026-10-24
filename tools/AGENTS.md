# Script organization

- Keep application and utility scripts directly in this directory; do not create application subdirectories.
- Keep all ADNI-specific function definitions in `adni-tools.R`; ADNI entrypoints should only load it and invoke the requested action.
- Keep all Ocean-specific function definitions in `ocean-tools.R`; Ocean entrypoints should only load it and invoke the requested action.
- Keep reusable functions shared by ADNI, Ocean, and other analyses in the package.
- Keep report Rmd/HTML files and manuscript assets under `paper/`, and data under `data/`.
