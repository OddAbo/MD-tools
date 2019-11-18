! 2019-11-14 by Dianmo
!
! Radial distribution function g(r)

program rdf
  implicit none
  real, parameter :: dt = 0.002, box = 20, pi = 3.1415926
  integer, parameter :: natom = 43, ngr = 300
  logical :: iexist
  integer :: nstep, istep, igr, tstart, tend, i, j, k
  real(kind=8) :: rv(6,natom), rho, l_box
  real(kind=8) :: dgr, dr(3), dist, rcut2, g(ngr), vgr

  l_box = 10.
  rho = float(natom) / l_box**3
  dgr =  0.5 * box / ngr
  g = 0
  rcut2 = (box/2) ** 2

  inquire(file="md.out",exist=iexist)
  if(.not.iexist) then
    stop "Error: ""md.out"" doesn't exist!"
  endif

  inquire(file="rdf.dat",exist=iexist)
  if(iexist) then
    stop "Error: ""rdf.dat"" already exists!&
    & Rename or remove it before running this program again."
  endif
  
  open(10,file="md.out",form="unformatted",status="old")
  open(20,file="rdf.dat",form="formatted",status="new")
  rewind(10)
  rewind(20)

  write(*,*) "Starts at time(ps):"
  read(*,*) tstart
  write(*,*) "Ends at time(ps):"
  read(*,*) tend

  if (tstart > tend) then
    write(*,*) "Error: End time must be later than start time!"
    stop
  endif
  
  nstep = int((tend-tstart)/dt)

  skip: do i = 2, floor(tstart/dt)
    read(10)
  enddo skip

  do i = 1, nstep
    read(10) istep, rv

    atom1: do j = 1, natom-1
      atom2: do k = j+1, natom
        dr = rv(1:3,j) - rv(1:3,k)

        ! periodic boundary condition
        ! dr = dr - box * nint(dr/box)

        dist = sum(dr**2)
        if (dist > rcut2) then
          cycle atom2
        endif
        dist = sqrt(dist)
        igr = int(dist/dgr)
        g(igr) = g(igr) + 2

      enddo atom2
    enddo atom1
  enddo

  ! Attention !
  ! To cluster, density-based correction couldn't be made 
  ! due its irregular shape.
  do i = 1, ngr
    vgr = (4 * pi) / 3. * ((i+1)**3 - i**3) * dgr**3
    g(i) = g(i) / (natom * nstep * vgr * rho)
    write(20,"(2(f10.6))") (i-0.5)*dgr, g(i)
  enddo

end program rdf