# Little tools for MD analysing

*To avoid overwriting the original file, almost all program would require to rename or remove it before running again.*

## [diffusion.f90](https://github.com/OddAbo/MD-tools/blob/master/diffusion.f90)
`plot 'diff.dat' u 1:2 w l`
\t  mean-squared displacement

`plot 'diff.dat' u 3:4 w l`
\t  velocity auto-correlation function

`plot 'diff.dat' u 5:6 w l`
\t  phonon density of states