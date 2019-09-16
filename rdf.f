      ! Radial distribution function g(r) of particles
      ! 2019-07-11
      ! Author: dianmo

      program rdf

        implicit none
        
        real, parameter :: Rcutoff = 10
        real, parameter :: Xlength = 32.72, Ylength = 32.72, Zlength = 32.72
        integer, parameter :: Natom = 43, Maxatom = 43, Nbox = 300
        integer :: istep, ipbc, nstep, i, j, k, l, ibox
        real*8 :: rv(6,Maxatom), dr(3), volume
        real*8 :: distance, g(Nbox), Rcut2, dbox
        
        nstep = 10000                      ! takes n steps into account
        ipbc = 0                           ! periodic boundary condition
        dbox = Rcutoff / Nbox              ! divided Rcut into several boxes
        g = 0.                             ! initializes g(r)
        Rcut2 = Rcutoff**2
  
        
        open(01,file='md.out',form='unformatted',status='unknown')
        open(02,file='rdf.dat',form='formatted',status='unknown')
        rewind(01)
        rewind(02)
        
        do l = 1, 5
  
          skip: do i = 1, 30000
            read(01)
          enddo skip
          
          do i = 1, nstep
            read(01) istep, rv
  
            part1: do j = 1, Natom-1
              part2: do k = j+1, Natom
                
                dr = (rv(1:3,j) - rv(1:3,k))
                
                pbc: if (ipbc == 1) then
                      dr(1) = Xlength * (dr(1) - nint(dr(1)))
                      dr(2) = Ylength * (dr(2) - nint(dr(2)))
                      dr(3) = Zlength * (dr(3) - nint(dr(3)))
                end if pbc
  
                distance = sum(dr**2)
                if (distance - Rcut2 > 0) then
                  cycle part2
                else
                  distance = sqrt(distance)
                  ibox = (distance / dbox)
                  g(ibox) = g(ibox) + 2
  
              enddo part2
            enddo part1
          enddo
        enddo
          
        do i = 1, Nbox-1
          volume = ((i+1)**3-i**3)*dbox**3
          g(i) = g(i) / (Natom*nstep*volume) / 5.     ! normalization
          write(02,'(2(2XF10.6))') (i+0.5)*dbox, g(i)
        enddo
  
        end program rdf
  