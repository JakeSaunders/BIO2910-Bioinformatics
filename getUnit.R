## getUnit() -- download one unit folder from the BIO 2910 GitHub repo.
## Base R only: nothing to install, no git, no GitHub account. Public repos only.
##   getUnit("Unit02","BIO2910")   ->  ./BIO2910/Unit02

getUnit <- function( unit, # folder to pull from the repo, e.g. "Unit02"
                    dest   = "BIO2910",   # folder to put the unit into
                    repo   = "JakeSaunders/BIO2910-Bioinformatics",
                    branch = "main") {

  # two unused temp paths: one for the zip, one to unpack into
  z <- tempfile(fileext = ".zip"); d <- tempfile()

  # every public repo is downloadable at this URL; mode="wb" or the zip corrupts
  download.file(sprintf("https://github.com/%s/archive/refs/heads/%s.zip", repo, branch),
                z, mode = "wb")

  unzip(z, exdir = d)

  # GitHub wraps everything in one "<repo>-<branch>" folder, so grab that, then the unit
  src <- file.path(list.dirs(d, recursive = FALSE)[1], unit)

  # catch typos here; file.copy() would just return FALSE with a vague warning
  if (!dir.exists(src)) stop("'", unit, "' is not a folder in this repo.")

  # file.copy(recursive=TRUE) needs dest to exist already; it won't create it
  if (!dir.exists(dest)) dir.create(dest, recursive = TRUE)

  # copies the folder INTO dest (result is dest/unit), with all subfolders
  file.copy(src, dest, recursive = TRUE)

  unlink(c(z, d), recursive = TRUE)       # delete the temp zip and unpacked copy

  out <- normalizePath(file.path(dest, unit))   # relative path -> full path
  message("Downloaded ", unit, " -> ", out)
  invisible(out)                          # returns the path without printing it
}
