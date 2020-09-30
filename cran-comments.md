## Comments for CRAN submission

### Release
This is a fourth minor release (0.4.0).

### Test environments
* Windows 10 version 2004, R 4.0.2 (local)
* Ubuntu Linux 16.04 LTS, R 4.0.2 (Travis CI)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (R-hub builder)
* Ubuntu Linux 16.04 LTS, R-release, GCC (R-hub builder)
* Fedora Linux, R-devel, clang, gfortran (R-hub builder)

### R CMD check results

#### NOTES:

1. checking CRAN incoming feasibility ... NOTE
Maintainer: 'Mark Murphy <murphymv@gmail.com>'
Found the following (possibly) invalid URLs:
  URL: https://doi.org/b8b782
    From: man/R2.Rd
    Status: 403
    Message: Forbidden
  URL: https://doi.org/bvxb6s
    From: man/glt.Rd
    Status: Error
    Message: libcurl error code 56:
      	Send failure: Connection was reset

   \- These URLs open fine when tested locally (multiple browsers).

2. checking installed package size ... NOTE
  installed size is  7.0Mb
  sub-directories of 1Mb or more:
    data   6.7Mb

   \- Example data provided with this package is used to facilitate the running 
   of quick examples, and every effort has been made to minimise file size 
   without loss of demonstrative value to the user.

3. checking for future file timestamps ... NOTE
unable to verify current time
   
   \- Apparently an issue with worldclockapi.com: 
   https://stackoverflow.com/questions/63613301/r-cmd-check-note-unable-to-verify-current-time.

4. checking examples ... NOTE
Examples with CPU or elapsed time > 5s
       user system elapsed
semEff 6.36  0.077   6.438

   \- Testing of examples apparently only exceeds recommended time limits on 
   Linux systems (R-hub builder).

### Downstream dependencies
There are currently no downstream dependencies.
