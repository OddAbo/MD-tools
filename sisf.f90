! 2019-11-21 by Dianmo
!
! The space Fourier transform of the self part of the van Hove correlation
! function G_s(r,t) gives the self-intermediate-scattering function F_s(q,t)
! (also called auto-correlation function).

program sisf
  implicit none
  logical :: iexist
  real, parameter :: pi = 3.14159265, dt = 0.002
  real, parameter :: rmax = 10, dr = 0.1
  integer, parameter :: natom = 43, nstep = 1000000
  integer :: iend, iatom1, iatom2, istep, ibox, nbox, i, j
  real(kind=8) :: rv_0(6,natom), rv_new(6,natom), qmax, rdamp
  real(kind=8) :: fs, dist, rmax2, gs(100)    ! gs(nbox)
  
  nbox = int(rmax / dr)
  rmax2 = rmax * rmax
  qmax = 2 * pi / 2.7
  rdamp = rmax / sqrt(3.0)
  inquire(file="md.out",exist=iexist)
  if(.not.iexist) then
    stop "Error: ""md.out"" doesn't exist!"
  endif

  ! inquire(file="vhcf.dat",exist=iexist)
  ! if(iexist) then
  !   stop "Error: ""vhcf.dat"" already exists!&
  !   & Rename or remove it before running this program again."
  !   stop
  ! endif

  open(10,file="md.out",form="unformatted",status="old")
  open(20,file="vhcf.dat", form="formatted",status="unknown")
  rewind(10)
  rewind(20)
  write(20,*) "t      r      Gs(r,t)"
  
  read(10) istep, rv_0
  rewind(10)

  loop: do i = 1, nstep
    fs = 0
    read(10,iostat=iend) istep, rv_new
    if(is_iostat_end(iend)) exit loop

    ! Self part of van Hove correlation function
    gs = 0
    particle1: do iatom1 = 1, natom
      particle2: do iatom2 = 1, natom
        dist = sum((rv_new(1:3,iatom1) - rv_0(1:3,iatom1))**2)
        ! if (dist > rmax2) cycle particle
        dist = sqrt(dist)
        ibox = int(dist/dr) + 1
        gs(ibox) = gs(ibox) + 1
      enddo particle2
    enddo particle1

    gs = gs / float(natom)

    ! Fs from space Fourier transform of Gs
    do j = 1, nbox
      fs = fs + gs(j) * cos((j-1) * dr * qmax) * dr! * exp(-((j-1)*dr/rdamp)**2)
    enddo
    write(20,"(2f10.6)") (i-1)*dt, fs
    if (fs < 0) exit loop
  enddo loop

end program sisf
