# Little tools for MD analysing

**To avoid overwriting the original file, almost all program would require to rename or remove it before running again.**

**You might need to make a little change to original code, like the initialization of array** ***rv***
```
real(kind=8) :: rv(6,natom)
```
## [diffusion.f90](https://github.com/OddAbo/MD-tools/blob/master/diffusion.f90)
```
plot 'diff.dat' u 1:2 w l        # mean-squared displacement
plot 'diff.dat' u 3:4 w l        # velocity auto-correlation function
plot 'diff.dat' u 5:6 w l        # phonon density of states
```