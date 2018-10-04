# FYS4150 Project 2 Unittests
Unittests.f08 contains two subroutines, one for testing of the module subroutine maxoffdiag and one for testing the module subroutine rotate. 

## Dependencies
gFortran

## Runs
Compilation:

$ gfortran -O3 -fopt-info -o unittests.exe unittests.f08 maxoffdiag_rotate.o

.exe:

C:\location_of_file>unittests.exe

Should return:

Subroutine maxoffdiag: All tests passed

Subroutine rotate: All tests passed


