! 2019-11-14 by Dianmo
!
! Structural overlap:
! Q(t) = θ(a-|ri(t)-rj(0))/N
! where θ(x) is the Heaviside function, a is cutoff radius
! Q1: i = j
! Q2: i=1,...,natom, j=1,...,natom

program overlap
  implicit none
  integer, parameter :: natom = 43, nblock = 100
  real, parameter :: dt = 0.002
  logical :: iexist
  integer :: istep, iread, i, j, q1, q2
  real(kind=8) :: rv(6,natom), pos_ref(3,natom), pos_new(3,natom)
  real(kind=8) :: dist2(natom,natom), r1, rcut2, t

  inquire(file="md.out",exist=iexist)
  if(.not.iexist) then
    write(*,*) "Error: ""md.out"" doesn't exist!"
    stop
  endif

  inquire(file="overlap.dat",exist=iexist)
  if(iexist) then
    write(*,*) "Error: ""overlap.dat"" already exists!"
    stop
  endif
  
  open(10,file="md.out",form="unformatted",status="old")
  open(20,file="overlap.dat",form="formatted",status="new")
  rewind(10)
  rewind(20)

  write(*,*) "Time of reference structure: (ps)"
  read(*,*) t
  
  ! r1: radius of the first peak in rdf
  r1 = 2.7
  iread = 0
  istep = int(t / dt)
  pos_ref = 0
  rcut2 = (0.333333 * r1) ** 2
  
  skip: do i = 1, istep
    read(10)
  enddo skip

  ref_struc: do i = 1, nblock
    read(10) istep, rv
    pos_ref = pos_ref + rv(1:3,1:natom)
  enddo ref_struc
  pos_ref = pos_ref / float(nblock)
  rewind(10)
  t = 0

  do while (.true.)
    pos_new = 0
    q1 = 0
    q2 = 0

    new_struc: do i = 1, nblock
      read(10,iostat=iread) istep, rv
      if (iread<0) stop
      pos_new = pos_new + rv(1:3,1:natom)
    enddo new_struc
    pos_new = pos_new / float(nblock)
    t = t + nblock

    distance: do i = 1, natom
      do j = 1, natom
        dist2(i,j) = sum((pos_ref(1:3,i) - pos_new(1:3,j)) ** 2)
      enddo
    enddo distance

    all_atom: do i = 1, natom
      overlap1: if (dist2(i,i) < rcut2) then
        q1 = q1 + 1
        cycle all_atom
      endif overlap1

      overlap2: do j = 1, natom
        if (i == j) cycle overlap2
        if (dist2(i,j) < rcut2) then
          q2 = q2 + 1
          exit overlap2
        endif
      enddo overlap2
    enddo all_atom

    write(20,"(3(f10.3))") t*dt, q1/float(natom), q2/float(natom)

  enddo

end program overlap