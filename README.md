# Little tools for MD analysing

## diffusion.f90:
> mean-squared displacement:

`plot 'diff.dat' u 1:2 w l`

> velocity auto-correlation function:

`plot 'diff.dat' u 3:4 w l`

> phonon density of states:

`plot 'diff.dat' u 5:6 w l`