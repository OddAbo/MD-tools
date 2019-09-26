      ! 2019-07-08
      ! Author: Dianmo Zhang

      ! Structural overlap: 
      ! Q(t) = θ(a-|r_i(t)-r_j(0)|)/N
      ! where θ(x) is the Heaviside function, a is cutoff radius
      ! Q(1): compares with their own initial position only
      ! Q(2): eliminates the effects of exchange motion

      program overlap
      
      implicit none

      integer, parameter :: Natom = 43, Maxatom = 43, Naverage = 500
      integer :: istep, i, j, l, Q(2), time
      real*8 :: rv(6,Maxatom), pos_0(3,Natom), pos_new(3,Natom)
      real*8 :: dr2(Natom,Natom), r_cutoff2, temperature

      time = 0
      r_cutoff2 = 0.685584
      
      open(01,file='md.out',form='unformatted',status='unknown')
      open(02,file='overlap.dat',form='formatted',status='unknown')
      rewind(01)
      rewind(02)

      ! jumps to desired moment
      skip: do i = 1, 25000
        read(01)
      enddo skip
      time = i

      ! stores coordinates of particles of inherent structure in pos_0
      inherent_struc: do i = 1, Naverage
        read(01) istep, rv
        pos_0 = pos_0 + rv(1:3,1:Natom)
      enddo inherent_struc
      pos_0 = pos_0 / float(Naverage)
      time = time + Naverage

      do while (.true.)
        Q = 0
        pos_new = 0.
        
        ! stores coordinates of particles of new structure in pos_new
        new_struc: do l = 1, Naverage
          read(01,end=999) istep, rv
          pos_new = pos_new + rv(1:3,1:Natom)
        enddo new_struc
        pos_new = pos_new / float(Naverage)
        time = time + Naverage
        
        distance: do i = 1, Natom
          do j = 1, Natom
            dr2(i,j) = sum((pos_new(1:3,i) - pos_0(1:3,j))**2)
          enddo
        enddo distance
        
        all_particle: do i = 1, Natom
          
          Q1: if (dr2(i,i) < r_cutoff2) then 
            Q = Q + 1
            cycle all_particle
          end if Q1

          Q2: do j = 1, Natom
            if (dr2(i,j) < r_cutoff2) then
              Q(2) = Q(2) + 1
              exit Q2
            end if
          enddo Q2

        enddo all_particle
          
        write(02,*) time*0.001, Q/float(Natom)
    
      enddo
999   continue

      end program overlap
