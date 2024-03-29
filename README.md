# Little tools for MD analyzing
> To avoid overwriting the original file due to misoperation, almost all program would require to rename or remove it before running again.

You might need to make a little change to original code, like **the shape of array** ***rv***
```
real(kind=8) :: rv(6,natom)
```

---
## [diffusion.f90](https://github.com/OddAbo/MD-tools/blob/master/diffusion.f90)
```
plot 'diff.dat' u 1:2 w l        # mean-squared displacement
plot 'diff.dat' u 3:4 w l        # velocity auto-correlation function
plot 'diff.dat' u 5:6 w l        # phonon density of states
```

---
## [eu-dist.f90](https://github.com/OddAbo/MD-tools/blob/master/eu-dist.f90)
```
plot 'eud.dat' u 1:2 w l        # Euclidean distance
```

---
## [overlap.f90](https://github.com/OddAbo/MD-tools/blob/master/overlap.f90)
```
plot 'overlap.dat' u 1:2 w l        # Q1
plot 'overlap.dat' u 1:3 w l        # Q2
```

---
## [rdf.f90](https://github.com/OddAbo/MD-tools/blob/master/rdf.f90)

You might need to make a little change to **PBC** part if necessary
```
plot 'rdf.dat' u 1:2 w l        # radial distribution function
```

---
## [transform.f90](https://github.com/OddAbo/MD-tools/blob/master/transform.f90)

Transform specific fragment of data in the unformatted file **"md.out"** into formatted file **"md.dat"**
```
Starts at time(ps):
20
Ends at time(ps):
30
```

---
## [sisf.f90](https://github.com/OddAbo/MD-tools/blob/master/sisf.f90)

**PBC** hasn't been added into code
```
set log x 10
plot 'sisf.dat' u 1:2 w l        # self-intermediate scattering function
```

---
To be continue...