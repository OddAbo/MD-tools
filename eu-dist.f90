! 2019-11-13 by Dianmo
! 
! Euclidean distance

program eudist
  implicit none
  logical :: iexist
  integer, parameter :: natom = 43
  real, parameter :: dt = 0.002
  integer :: iend, istep, i, j, t
  real(kind=8) :: rv(6,natom), r(3,52500,natom)
  real(kind=8) :: dr2(natom), dr_bar

  t = 2500

  inquire(file="md.out",exist=iexist)
  if(.not.iexist) then
    stop "Error: ""md.out"" doesn't exist!"
  endif

  inquire(file="eud.dat",exist=iexist)
  if(iexist) then
    stop "Error: ""eud.dat"" already exists!&
    & Rename or remove it before running this program again."
  endif
  
  open(10,file="md.out",form="unformatted",status="old")
  open(20,file="eud.dat",form="formatted",status="new")
  rewind(10)
  rewind(20)

  divide: do while (.true.)

    store: do i = 1, 52500
      read(10,iostat=iend) istep, rv
      if(is_iostat_end(iend)) stop
      r(1:3,i,1:natom) = rv(1:3,1:natom)
    enddo store

    cal: do i = 1, 50000
      dr_bar = 0

      all_atom: do j = 1, natom
        dr2(j) = sum((r(1:3,i,j) - r(1:3,i+2500,j))**2)
        dr_bar = dr_bar + dr2(j)
      enddo all_atom

      dr_bar = sqrt(dr_bar/float(natom))
      t = t + 1
      
      write(20,"(f8.3)",advance="no") t*dt
      ! outputs displacement of all atoms:
      ! do j = 1, natom
      !   write(20,"(f10.6)",advance="no") dr2(j)
      ! enddo
      write(20,"(f10.6)") dr_bar
    enddo cal

    back: do i = 1, 2500
      backspace (10)
    enddo back

  enddo divide

end program eudist