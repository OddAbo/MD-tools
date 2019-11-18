! 2019-11-09 by Dianmo
! 
! This small program calculates:
! mean-squared displacement,
! velocity auto-correlation funciton
! and phonon density of state.
!
! See more in ./reference/VAC.pdf

program diffusion
  implicit none
  integer, parameter :: natom = 43, nstep = 2000000
  real, parameter :: dt = 0.002, PI = 3.14159265, tau = 2.886751346
  logical :: iexist
  integer :: istep, iatom, iend, i, j
  real(kind=8) :: msd(2500), vac(2500), dos(2500), omega(2500)
  real(kind=8) :: t, d_omega
  real(kind=8) :: rv(6,natom), rv_new(6,natom)

  inquire(file="md.out",exist=iexist)
  if(.not.iexist) then
    stop "Error: ""md.out"" doesn't exist!"
  endif

  inquire(file="diff.dat",exist=iexist)
  if(iexist) then
    stop "Error: ""diff.dat"" already exists!&
    & Rename or remove it before running this program again."
    stop
  endif

  open(10,file="md.out",form="unformatted",status="old")
  open(20,file="diff.dat", form="formatted",status="new")
  rewind(10)
  rewind(20)
  write(20,*) "t    msd    t    vac    omega    dos"
  msd = 0
  vac = 0
  dos = 0
  omega = 0
  iend = 0
  d_omega = 20 * PI / 2500.

  do while (.true.)
    read(10,iostat=iend) istep, rv
    if(is_iostat_end(iend)) stop
    do i = 2, 2500
      read(10,iostat=iend) istep, rv_new
      if(is_iostat_end(iend)) stop
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
    write(20,"(3(f8.3,f12.6))") &
    i*0.002, msd(i), i*0.002, vac(i), i*d_omega/PI/2., omega(i)
  enddo

end program diffusion