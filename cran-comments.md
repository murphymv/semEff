## Release
This is a second version release (0.2.0)

## Test environments
* Windows 10 version 1909, R 3.6.2 (local)
* Ubuntu Linux 16.04 LTS, R 3.6.1 (travis-ci)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (R-hub builder)
* Ubuntu Linux 16.04 LTS, R-release, GCC (R-hub builder)
* Fedora Linux, R-devel, clang, gfortran (R-hub builder)

## R CMD check results
There were two NOTES:

* checking CRAN incoming feasibility ...NB: need Internet access to use CRAN incoming checks
 NOTE
Maintainer: ‘Mark Murphy <murphymv@gmail.com>’

Possibly mis-spelled words in DESCRIPTION:
  Lefcheck (10:63)

Found the following (possibly) invalid URLs:
  URL: https://doi.org/b8b782
    From: man/R2.Rd
    Status: 403
    Message: Forbidden
  URL: https://doi.org/bvxb6s
    From: man/getY.Rd
    Status: Error
    Message: libcurl error code 56:
      	OpenSSL SSL_read: SSL_ERROR_SYSCALL, errno 104

* checking installed package size ... NOTE
  installed size is  7.8Mb
  sub-directories of 1Mb or more:
    data   7.5Mb
    
The word in DESCRIPTION is an author's name. The URLs open fine when tested
locally (multiple browsers).

* checking installed package size ... NOTE
  installed size is  7.8Mb
  sub-directories of 1Mb or more:
    data   7.5Mb

Example data provided with this package is used to facilitate the running of
quick examples, and every effort has been made to minimise file size without
loss of demonstrative value to the end user.

## Downstream dependencies
There are currently no downstream dependencies for this package.
