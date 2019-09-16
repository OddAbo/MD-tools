      ! Saves data in 'md.out' to a formatted file
      ! 2019-09-10
      ! Author: dianmo

      program format

      implicit none
      
      integer, parameter :: natom = 13
      real*8 rv(6,natom)
      integer i

      open(01,file='md.out',form='unformatted',status='unknown')
      open(02,file='formatted.dat',form='formatted',status='unknown')
      rewind(01)
      rewind(02)

      do while (.true.)
        read(01,end=999) i, rv
        write(02,*) rv(1:3,:)
      enddo

999   continue

      end program format
