## Release
This is a fourth minor release (0.4.0).

## Test environments
* Windows 10 version 2004, R 4.0.2 (local)
* Ubuntu Linux 16.04 LTS, R 4.0.2 (travis-ci)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (R-hub builder)
* Ubuntu Linux 16.04 LTS, R-release, GCC (R-hub builder)
* Fedora Linux, R-devel, clang, gfortran (R-hub builder)

## R CMD check results
NOTES:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Mark Murphy <murphymv@gmail.com>’

Found the following (possibly) invalid URLs:
  URL: https://doi.org/b8b782
    From: man/R2.Rd
    Status: 403
    Message: Forbidden
  URL: https://doi.org/bvxb6s
    From: man/glt.Rd
    Status: Error
    Message: libcurl error code 56:
      	OpenSSL SSL_read: SSL_ERROR_SYSCALL, errno 104

These URLs open fine when tested locally (multiple browsers).

* checking installed package size ... NOTE
  installed size is  7.0Mb
  sub-directories of 1Mb or more:
    data   6.7Mb

Example data provided with this package is used to facilitate the running of
quick examples, and every effort has been made to minimise file size without
loss of demonstrative value to the end user.

* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
        user system elapsed
predEff 3.24  0.024   9.553

Testing of examples apparently only exceeds recommended time limits on Linux
systems (R-hub builder)?

## Downstream dependencies
There are currently no downstream dependencies.
