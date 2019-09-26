      ! diffusion coefficient 
      ! & velocity correlation funtion
      ! & density of frequency
      ! 2019-09-24
      ! Author: dianmo

      program diffusion
      
      implicit none

      integer, parameter :: npart = 13
      real, parameter :: PI = 3.14159265
      real :: msd(2500), vcf(2500), dos(2500), omega(2500)
      real :: t, tau, r2, v2, d_omega
      real*8 :: rv(6,npart), rv_new(6,npart)
      integer :: istep, ipart, ivcf, idos
      integer :: i, j, k

      open(01,file='md.out',form='unformatted',status='unknown')
      open(02,file='dos.dat',form='formatted',status='unknown')
      open(03,file='new.dat',form='formatted',status='unknown')
      rewind(01)
      rewind(02)
      rewind(03)

      ivcf = 1
      idos = 1
      msd = 0
      vcf = 0
      dos = 0
      tau = 2.886751346                    ! tau = t_max / sqrt(3)
      omega = 0
      d_omega = 20 * PI / 2500.
      
      whole_file: do while (.true.)
        read(01,end = 999) istep, rv
        
        block_average: do i = 2, 2500
          read(01,end = 999) istep, rv_new
          r2 = 0
          v2 = 0

          all_part: do ipart = 1, npart
            msd(i) = msd(i) + sum((rv(1:3,ipart)-rv_new(1:3,ipart))**2)
            if (ivcf == 1) vcf(i) = vcf(i) + 
     &      sum(rv(4:6,ipart)*rv_new(4:6,ipart)) / sum(rv(4:6,ipart)**2)
          enddo all_part
        enddo block_average
      enddo whole_file
999   continue

      ! normalization
      msd = msd / float(nstep * npart) * 2500.
      vcf = vcf / float(nstep * npart) * 2500.
      vcf(1) = 1

      if (idos == 1) then
        do i = 1, 2500
          do j = 1, 2500
            t = j * 0.002
            omega(i) = omega(i) + 
     &      vcf(j) * exp(-(t/tau)**2) * cos(i*d_omega*t) * 0.004
          enddo
        enddo
      end if
      
      do i = 1, 2500
        write(02,'(6(f12.8))') 
     &    i*0.002, msd(i), i*0.002, vcf(i), i*d_omega/PI/2., omega(i)
      enddo

      end program diffusion