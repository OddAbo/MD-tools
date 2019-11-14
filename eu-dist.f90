! 2019-11-13 by Dianmo
! 
! Euclidean distance

program eudist
  implicit none
  logical :: iexist
  integer, parameter :: natom = 43, nstep = 1000000
  real, parameter :: dt = 0.002
  integer :: istep, i, j, k, t
  real(kind=8) :: rv(6,natom), r(3,52500,natom)
  real(kind=4) :: dr2(natom), msd_tot

  t = 2500

  inquire(file="md.out",exist=iexist)
  if(.not.iexist) then
    write(*,*) "Error: ""md.out"" doesn't exist!"
    stop
  endif

  inquire(file="eud.dat",exist=iexist)
  if(iexist) then
    write(*,*) "Error: ""eud.dat"" already exists!"
    stop
  endif
  
  open(10,file="md.out",form="unformatted",status="old")
  open(20,file="eud.dat",form="formatted",status="new")
  rewind(10)
  rewind(20)

  divide: do i = 1, (nstep-2500)/50000

    store: do j = 1, 52500
      read(10) istep, rv
      r(1:3,j,1:natom) = rv(1:3,1:natom)
    enddo store

    cal: do j = 1, 50000
      msd_tot = 0

      all_atom: do k = 1, natom
        dr2(k) = sum((r(1:3,j,k) - r(1:3,j+2500,k))**2)
        msd_tot = msd_tot + dr2(k)
      enddo all_atom

      msd_tot = sqrt(msd_tot/float(natom))
      t = t + 1
      
      write(20,"(f12.3)",advance="no") t*dt
      ! outputs msd of all atoms:
      ! do k = 1, natom
      !   write(20,"(f12.6)",advance="no") dr2(k)
      ! enddo
      write(20,"(f12.6)") msd_tot
    enddo cal

    back: do j = 1, 2500
      backspace (10)
    enddo back
  enddo divide
end program eudist