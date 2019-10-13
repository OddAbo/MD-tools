      program format

      implicit none

      integer, parameter :: natom = 43
      real*8 rv(6,natom)
      integer i, j

      open(01,file='md.out',form='unformatted',status='unknown')
      open(02,file='position.dat',form='formatted',status='unknown')
      rewind(01)
      rewind(02)

      do while (.true.)
        read(01,end=999) i, rv
        write(02,'(129(F15.9))') rv(1:3,:)
      enddo
999   continue

      end program format
