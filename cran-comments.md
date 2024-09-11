## Comments for CRAN submission

### Release

This is a resubmission (0.7.1).

This version corrects the following issues:

NOTE

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10/bqd43d (moved to https://doi.org/10.1890/08-1034.1)
    From: README.md
    Status: 301
    Message: Moved Permanently
  URL: https://doi.org/10/cqm32d (moved to https://doi.org/10.1207/S15328007SEM0702_4)
    From: README.md
    Status: 301
    Message: Moved Permanently
  URL: https://doi.org/10/f8s8rb (moved to https://doi.org/10.1111/2041-210X.12512)
    From: README.md
    Status: 301
    Message: Moved Permanently
  URL: https://doi.org/bqd43d (moved to https://doi.org/10.1890/08-1034.1)
    From: inst/doc/predicting-effects.html
    Status: 301
    Message: Moved Permanently
  URL: https://doi.org/d8gvwm (moved to https://doi.org/10.1890/1051-0761(2006)016%5B0503:ASEMAO%5D2.0.CO;2)
    From: inst/doc/semEff.html
    Status: 301
    Message: Moved Permanently
  URL: https://doi.org/f8s8rb (moved to https://doi.org/10.1111/2041-210X.12512)
    From: inst/doc/semEff.html
    Status: 301
    Message: Moved Permanently
  URL: https://www.buymeacoffee.com/murphymv (moved to https://buymeacoffee.com/murphymv)
    From: DESCRIPTION
          man/semEff-package.Rd
          README.md
    Status: 301
    Message: Moved Permanently

FIX

The buymeacoffee.com link has been fixed. The DOI links are https://shortdoi.org/ shortcuts to the original DOIs and are all valid (checked) and have been accepted in previous submissions. Is it possible to add these to a whitelist?

NOTE

Found the following Rd file(s) with Rd \link{} targets missing package
anchors:
  R2.Rd: predict.merMod, hatvalues.merMod
  bootCI.Rd: boot.ci
  bootEff.Rd: bootMer, boot
  predEff.Rd: predict.merMod, boot.ci
  semEff.Rd: boot.ci

FIX

Package anchors have been added where appropriate.

### R CMD check results

0 errors | 0 warnings | 0 notes

### revdepcheck results

I checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

-   I saw 0 new problems
-   I failed to check 0 packages

### Test environments

| System                                    | Source                           | R version  |
|--------------------------------|-------------------------|----------------|
| Windows 10 Home 22H2                      | Local                            | R 4.4.1    |
| Windows Server 2022 10.0.20348 Datacenter | Remote (GitHub Actions workflow) | R-release  |
| Mac OS 14.6.1 23G93                       | Remote (GitHub Actions workflow) | R-release  |
| Ubuntu 22.04.4 LTS                        | Remote (GitHub Actions workflow) | R-devel    |
| Ubuntu 22.04.4 LTS                        | Remote (GitHub Actions workflow) | R-release  |
| Ubuntu 22.04.4 LTS                        | Remote (GitHub Actions workflow) | R-oldrel-1 |
