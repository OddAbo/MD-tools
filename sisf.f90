! 2019-11-21 by dianmo
! Self-intermediate-scattering function Fs(q,t)

program sisf
  implicit none
  logical alive
  real, parameter :: pi = 3.14159265, dt = 0.002, tskip = 50
  integer, parameter :: natom = 43, maxatom = 43, nstep = 5000000, nstart = 100
  integer :: istart, istep, iatom, length_line, length_rec, irec, i
  integer :: fcount(nstep), stepmax
  real(kind=8) :: rv_ref(6,maxatom), rv_new(6,maxatom), Fs(nstep)
  real(kind=8) :: qmax, dr, Fs_tmp

  ! 2.7 is the radius of the first peak in rdf
  qmax = 2 * pi / 2.7
  fcount = 0
  stepmax = 0
  Fs = 0
  length_rec = 3*sizeof(istep) + sizeof(Fs_tmp)*size(rv_ref)
  ! If you are using Intel's compiler "ifort", uncomment the code below:
  ! length_rec = length_rec / 4
  
  inquire(file="md.out",exist=alive)
  if (.not.alive) then
    stop "Error: ""md.out"" doesn't exist!"
  endif

  open(10,file="md.out",action="read",access="direct",recl=length_rec)
  open(20,file="sisf.dat",action="write",form="formatted",status="new")
  rewind(20)

  new_start: do istart = 1, nstart
    ! sets start structure as reference
    irec = (istart-1)*(tskip/dt) + 1
    write(*,*) irec
    read(10,rec=irec) length_line, istep, rv_ref, length_line
    i = 1

    cal_Fs: do while (irec < nstep)
      Fs_tmp = 0
      read(10,rec=irec) length_line, istep, rv_new, length_line

      all_atom: do iatom = 1, natom
        ! TODO: periodic boundary condition
        dr = sum((rv_new(1:3,iatom)-rv_ref(1:3,iatom))**2)
        dr = sqrt(dr)
        Fs_tmp = Fs_tmp + cos(qmax*dr)
      enddo all_atom
      
      Fs(i) = Fs(i) + Fs_tmp
      fcount(i) = fcount(i) + 1
      irec = irec + 1
      stepmax = max(stepmax,i)
      i = i + 1
      if (Fs_tmp < 0) exit cal_Fs
    enddo cal_Fs
  enddo new_start
  close(10)

  Fs(1:stepmax) = Fs(1:stepmax) / nstart / natom

  print: do i = 1, stepmax
    write(20,"(f10.3,f10.6)") (i-1)*dt, Fs(i)
    if (Fs(i) < 0.0001) exit print
  enddo print
  close(20)

end program sisf
