! 2019-11-09 by Dianmo
! 
! This small program calculates:
! diffusion coefficient,
! velocity auto-correlation funciton
! and phonon density of state.

program diffusion
  implicit none
  integer, parameter :: natom = 43, nstep = 2000000
  real, parameter :: dt = 0.002, PI = 3.14159265, tau = 2.886751346
  logical :: iexist
  integer :: istep, iatom, iread, i, j
  real :: msd(2500), vac(2500), dos(2500), omega(2500)
  real :: t, d_omega
  real(kind=8) :: rv(6,natom), rv_new(6,natom)

  inquire(file="md.out",exist=iexist)
  if(.not.iexist) then
    write(*,*) "Error: ""md.out"" doesn't exist!"
    stop
  endif

  inquire(file="diff.dat",exist=iexist)
  if(iexist) then
    write(*,*) "Error: ""diff.dat"" already exists!"
    stop
  endif

  open(10,file="md.out",form="unformatted",status="old")
  open(20,file="diff.dat", form="formatted",status="new")
  rewind(10)
  rewind(20)

  msd = 0
  vac = 0
  dos = 0
  omega = 0
  iread = 0
  d_omega = 20 * PI / 2500.

  do while (iread==0)
    read(10,iostat=iread) istep, rv

    do i = 2, 2500
      read(10,iostat=iread) istep, rv_new
      msd(i) = msd(i) + sum((rv(1:3,:) - rv_new(1:3,:))**2)

      do iatom = 1, natom
        vac(i) = vac(i) + &
        sum(rv(4:6,iatom)*rv_new(4:6,iatom)) / sum(rv(4:6,iatom)**2)
      enddo

    enddo

  enddo

  msd = msd / float(nstep * natom) * 2500.
  vac = vac / float(nstep * natom) * 2500.
  vac(1) = 1
  
  do i = 1, 2500
    do j = 1, 2500
      t = j * dt
      omega(i) = omega(i) + &
      0.004 * vac(j) * exp(-(t/tau)**2) * cos(i*d_omega*t)
    enddo
  enddo

  do i = 1, 2500
    write(20,"(6(F12.6))") &
    i*0.002, msd(i), i*0.002, vac(i), i*d_omega/PI/2., omega(i)
  enddo

end program diffusion