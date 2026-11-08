#!/usr/bin/env Rscript
# Resume the frozen OAI eivGP workflow using the current implementation.
if (sys.nframe() == 0L) {
  script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]))
  source(file.path(dirname(script), 'resume-oai-tools.R'))
  resume_oai_main('eivGP', script)
}
