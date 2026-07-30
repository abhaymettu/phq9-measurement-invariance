#!/usr/bin/env Rscript
# Fetch the two NHANES 2017-2018 files this project needs and verify the vintage.
#
# Trap worth naming: the older-style path
#   https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DPQ_J.XPT
# still returns HTTP 200, but the body is a 20,905-byte CDC "Page Not Found" HTML
# page, not an XPT. A pipeline that only checks the status code will happily hand
# an HTML file to read_xpt. The live path is under /Nchs/Data/Nhanes/Public/.
# We check the SAS transport magic bytes, not the status code.

base <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles"

files <- list(
  DPQ_J.xpt  = "1c075f332e79736bd2360b41b0c10ed526883417e0b5996c407b1c564b01ca24",
  DEMO_J.xpt = "c0b46e0345ea19404928656277c8b0d10b0cca348a9b2fe4fc3c67e8b7ee73ec"
)

dir.create("data", showWarnings = FALSE)

for (f in names(files)) {
  dest <- file.path("data", f)
  if (!file.exists(dest)) {
    url <- file.path(base, f)
    cat("downloading", url, "\n")
    if (download.file(url, dest, quiet = TRUE) != 0) stop("download failed: ", url)
  }

  # SAS transport files start with this ASCII header record.
  magic <- rawToChar(readBin(dest, "raw", n = 26))
  if (!grepl("^HEADER RECORD", magic)) {
    stop(dest, " is not a SAS transport file. First bytes: ", substr(magic, 1, 26),
         "\nThe CDC serves its 404 page with status 200; check the URL path.")
  }

  got <- tools::sha256sum(dest)[[1]]
  cat(sprintf("%-12s %8d bytes  sha256 %s  %s\n", f, file.size(dest), substr(got, 1, 16),
              if (identical(got, files[[f]])) "vintage ok" else "SHA MISMATCH"))
  if (!identical(got, files[[f]]))
    warning("sha256 differs from the pinned vintage; CDC may have reissued the file.")
}
