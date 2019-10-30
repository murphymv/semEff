## Resubmission
This is a resubmission. The following issues were fixed:

* Removed quotation marks from the word 'piecewise' in the DESCRIPTION file.

* Added a reference for 'piecewise' structural equation models in the
description field of the DESCRIPTION file.

* Replaced \dontrun{} with \donttest{} in all Rd-files.

* Set no. of cores to two for all examples using parallel processing ("ncpus =
2").

## Release
This is a first submission to CRAN.

## Test environments
* Windows 10 version 1903, R 3.6.1 (local)
* Ubuntu Linux 16.04 LTS, R 3.6.1 (travis-ci)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (R-hub builder)
* Ubuntu Linux 16.04 LTS, R-release, GCC (R-hub builder)
* Fedora Linux, R-devel, clang, gfortran (R-hub builder)

## R CMD check results
There were two NOTES:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Mark Murphy <murphymv@gmail.com>'

New submission

Found the following (possibly) invalid URLs:
  URL: https://doi.org/b8b782
    From: man/R2.Rd
    Status: 403
    Message: Forbidden
  URL: https://doi.org/bvxb6s
    From: man/getY.Rd
    Status: Error
    Message: libcurl error code 56:
      	OpenSSL SSL_read: SSL_ERROR_SYSCALL, errno 10054

These links opened fine when tested locally (on multiple browsers).

* checking installed package size ... NOTE
  installed size is  7.8Mb
  sub-directories of 1Mb or more:
    data   7.5Mb

Example data provided with this package is used to facilitate the running of
quick examples, and every effort has been made to minimise file size without
loss of demonstrative value to the end user.

## Downstream dependencies
There are currently no downstream dependencies for this package.
