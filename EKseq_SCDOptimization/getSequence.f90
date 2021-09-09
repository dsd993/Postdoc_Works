module parammod

  implicit none

  integer :: num_res, num_seq, seq_length, rngseed
  character(len=1), dimension(:), allocatable :: resname
  integer, dimension(:), allocatable :: res_num
  real(8), dimension(:,:), allocatable :: hpres
  real(8), dimension(:), allocatable :: qres

end module parammod

program main

  use parammod
  use mtmod

  implicit none

  integer, parameter :: SEQLEN = 60
  integer, parameter :: NMAX = 1000000
  integer, parameter :: iseed = 3578921
  integer :: narg, nseq, ierr, num, nt
  character(len=255), dimension(:), allocatable :: args
  character(len=255) :: infile, line, hpfile, qfile
  character(len=1000) :: seq, seqo, seq1, seq2
  real(8), dimension(2,2) :: hp
  real(8) :: emin, emax, dx, x, e0, e0o, et, min_SCD, max_SCD
  character(len=1) :: ss
  integer :: imc, i, j, i1, i2
  logical :: arithmetic

  !below lines opens the 'input.in' file
  if (iargc() < 1) stop 'Need an input file!'
  call getarg(1, infile)
  open(10, file=trim(infile), status='old', iostat=ierr)
  if (ierr /= 0) stop 'Cannot open the input file!'
  num_res = 0
  num_seq = 0
  hpfile = ''
  qfile = ''
  arithmetic = .TRUE.
  
  !below loop read 'input.in' file to get the parameters for the sequences on wants to generate. 
  !here, for SCD optimization, pass 'q.dat' file that is read with the tag name 'qfile' and the file has info on residues and its charges.
  !'hp.dat' file that is read under the tag name 'hpfile' is required only for SHD calculation' 
  do while (ierr == 0)
     read(10, '(a)', iostat=ierr) line
     if (ierr == 0 .and. len_trim(line) > 0 .and. line(1:1) /= '#') then
        narg = LineArgs(line)
        allocate(args(narg))
        call GetLineArgs(line, narg, args)
        select case (trim(args(1)))
        case ('resname')
           num_res = narg - 1
           if (num_res < 1) stop 'No residue names!'
           allocate(resname(num_res))
           do i = 1, num_res
              resname(i) = trim(args(i+1))
           end do
        case ('resnum')
           if (num_res < 1) stop 'No residue names!'
           allocate(res_num(num_res))
           do i = 1, num_res
              read(args(i+1), *, iostat=ierr) res_num(i)
              if (ierr /= 0) stop 'Invalid format in resnum!'
           end do
        case ('numseq')
           read(args(2), *, iostat=ierr) num_seq
           if (ierr /= 0) stop 'Invalid format in numseq!'
        case ('seed')
           read(args(2), *, iostat=ierr) rngseed
           if (ierr /= 0) stop 'Invalid format in seed!'
        case ('hpfile')
           hpfile = trim(args(2))
           if (narg > 2) then
              if (trim(args(3)) == 'geometric') then
                 arithmetic = .FALSE.
              end if
           end if
        case ('qfile') 
           qfile = trim(args(2))
		case ('minSCD')
			read(args(2), *, iostat=ierr) min_SCD
			if (ierr /= 0) stop 'Invalid format in minSCD!'
		case ('maxSCD')
			read(args(2), *, iostat=ierr) max_SCD
			if (ierr /= 0) stop 'Invalid format in maxSCD!'
        end select
        deallocate(args)
     end if
  end do
  close(10)
 
  if (num_seq < 1) stop 'No sequence selected!'
  seq_length = sum(res_num)

  !allocating space for 1-D array 'qres' based on the number of residues 'num_res' in a seqeuence
  allocate(hpres(num_res,num_res))
  allocate(qres(num_res))
  hpres = 0.0
  qres = 0.0
  if (len_trim(hpfile) > 0) call getHPvalue(hpfile, arithmetic)
  if (len_trim(qfile) > 0) call getQvalue(qfile)

  call sgrnd(rngseed)

  !loop below forms the initial sequence to begin the optimization process
  !for eg: for E-K residues of 'num_res = 5' each, the initial sequence picked would be EEEEEKKKKK
  seq = ''
  nt = 0
  do i = 1, num_res
     do j = 1, res_num(i)
        nt = nt + 1
        seq(nt:nt) = resname(i)
     end do
  end do
  e0 = getSCD(seq)
  write(0,*) e0, trim(seq)
  
  if (min_SCD < 10000) then
    !below loop runs until the SCD value is pushed to the preferred minimum SCD value by the user
	do while (e0 < min_SCD)
		i1 = int(seq_length*grnd()) + 1
		i2 = int(seq_length*grnd()) + 1
		ss = seq(i1:i1)
		seq(i1:i1) = seq(i2:i2)
		seq(i2:i2) = ss
		e0 = getSCD(seq)
		write(0,*) e0, trim(seq)
	end do
	seq1 = seq
	emin = getSCD(seq1)
	write(0,*) "min=",emin, " ", trim(seq1)
  else
    !below loop runs for 1,000,000 iterations to pick the sequence with minimum SCD if there is no preference for minimum SCD
	do imc = 1, NMAX
		e0o = e0
		seqo = seq
		i1 = int(seq_length*grnd()) + 1
		i2 = int(seq_length*grnd()) + 1
		ss = seq(i1:i1)
		seq(i1:i1) = seq(i2:i2)
		seq(i2:i2) = ss
		e0 = getSCD(seq)
		if (e0 > e0o) then
			seq = seqo
			e0 = e0o
		end if
	end do
	seq1 = seq
	emin = getSCD(seq1)
	write(0,*) "min=",emin, " ", trim(seq1)
  end if
  
  !loop below forms the initial sequence to begin the optimization process
  !for eg: for E-K residues of 'num_res = 5' each, the initial sequence picked would be EEEEEKKKKK
  seq = ''
  nt = 0
  do i = 1, num_res
     do j = 1, res_num(i)
        nt = nt + 1
        seq(nt:nt) = resname(i)
     end do
  end do
  e0 = getSCD(seq)
  write(0,*) e0, trim(seq)
  
  if (max_SCD < 10000) then
    !below loop runs until the SCD value is pushed to the preferred maximum SCD value by the user
	do while (e0 < max_SCD)
		i1 = int(seq_length*grnd()) + 1
		i2 = int(seq_length*grnd()) + 1
		ss = seq(i1:i1)
		seq(i1:i1) = seq(i2:i2)
		seq(i2:i2) = ss
		e0 = getSCD(seq)
		write(0,*) e0, trim(seq)
	end do
	seq2 = seq
	emax = getSCD(seq2)
	write(0,*) "max=",emax, " ", trim(seq2)
  else
	!below loop runs for 1,000,000 iterations to pick the sequence with maximum SCD if there is no preference for maximum SCD by the user
    do imc = 1, NMAX
		e0o = e0
		seqo = seq
		i1 = int(seq_length*grnd()) + 1
		i2 = int(seq_length*grnd()) + 1
		ss = seq(i1:i1)
		seq(i1:i1) = seq(i2:i2)
		seq(i2:i2) = ss
		e0 = getSCD(seq)
		if (e0 < e0o) then
			seq = seqo
			e0 = e0o
		end if
	end do
	seq2 = seq
	emax = getSCD(seq2)
	write(0,*) "max=",emax, " ", trim(seq2)
  end if

  dx = (emax - emin)/dble(num_seq)
  write(0,*) emin, emax

  num = 1
  x = emin
  seq = seq1
  write(0,'(i8," ",a,f12.6)') num, trim(seq), getSCD(seq)
  write(6,'(i8," ",a,f12.6)') num, trim(seq), getSCD(seq)
  
  !below loop picks sequences with SCD values between the minimum and the maximum that is already picked
  do while (num < num_seq-1)
     x = x + dx
     et = getSCD(seq)
     e0 = (et - x)**2
     write(0,*) e0
     do imc = 1, NMAX
        e0o = e0
        seqo = seq
        i1 = int(seq_length*grnd()) + 1
        i2 = int(seq_length*grnd()) + 1
        ss = seq(i1:i1)
        seq(i1:i1) = seq(i2:i2)
        seq(i2:i2) = ss
        et = getSCD(seq)
        e0 = (et - x)**2
        if (e0 > e0o) then
           seq = seqo
           e0 = e0o
        end if
     end do
     write(0,*) e0
     num = num + 1
     write(0,'(i8," ",a,f12.6)') num, trim(seq), getSCD(seq)
     write(6,'(i8," ",a,f12.6)') num, trim(seq), getSCD(seq)
  end do
  num = num + 1
  seq = seq2
  write(0,'(i8," ",a,f12.6)') num, trim(seq), getSCD(seq)
  write(6,'(i8," ",a,f12.6)') num, trim(seq), getSCD(seq)
 
contains
  !below subroutine is only of interest for hydropathy calculation to calculate the lambda values
  subroutine getHPvalue(filename, comb)

    implicit none

    character(len=*), intent(in) :: filename
    logical, intent(in) :: comb
    integer :: ierr, narg, i, j
    real(8), dimension(num_res) :: hp_single
    character(len=255) :: line
    real(8) :: xx
    character(len=1) :: aa1, aa2
    logical :: single

    single = .FALSE.
    open(10, file=trim(filename), status='old', iostat=ierr)
    if (ierr /= 0) stop 'Cannot open the HP file!'
    do while (ierr == 0)
       read(10, '(a)', iostat=ierr) line
       if (ierr == 0 .and. len_trim(line) > 0) then
          narg = LineArgs(line)
          if (narg < 2) stop 'Invalid format in the HP file!'
          if (narg > 2) then
             read(line, *, iostat=ierr) aa1, aa2, xx
             if (ierr /= 0) stop 'Invalid format in the HP file!'
             i = SeqID(aa1)
             j = SeqID(aa2)
             if (i == 0 .or. j == 0) stop 'Missing residue!'
             hpres(i,j) = xx
             if (i /= j) hpres(j,i) = xx
          else
             read(line, *, iostat=ierr) aa1, xx
             if (ierr /= 0) stop 'Invalid format in the HP file!'
             i = SeqID(aa1)
             if (i == 0) stop 'Missing residue!'
             hp_single(i) = xx
          end if
       end if
    end do
    close(10)

    if (single) then
       do i = 1, num_res
          do j = 1, num_res
             if (comb) then
                hpres(i,j) = (hp_single(i) + hp_single(j))/2.0
             else
                hpres(i,j) = sqrt(hp_single(i)*hp_single(j))
             end if
          end do
       end do
    end if  

  end subroutine getHPvalue

  !below subroutine reads 'q.dat' file that contains the residues and its corresponding charges
  subroutine getQvalue(filename)

    implicit none

    character(len=*), intent(in) :: filename
    integer :: ierr, narg, i
    character(len=255) :: line
    real(8) :: xx
    character(len=1) :: aa1

    open(10, file=trim(filename), status='old', iostat=ierr)
    if (ierr /= 0) stop 'Cannot open the Q file!'
    do while (ierr == 0)
       read(10, '(a)', iostat=ierr) line
       if (ierr == 0 .and. len_trim(line) > 0) then
          narg = LineArgs(line)
          if (narg < 2) stop 'Invalid format in the Q file!'
          read(line, *, iostat=ierr) aa1, xx
          if (ierr /= 0) stop 'Invalid format in the Q file!'
          i = SeqID(aa1)
          if (i == 0) stop 'Missing residue!'
          qres(i) = xx
       end if
    end do
    close(10)

  end subroutine getQvalue

  !below function assigns sequence ids to the residues of interest
  integer function SeqID(s)

    implicit none

    character(len=1), intent(in) :: s
    integer :: i

    SeqID = 0
    do i = 1, num_res
       if (s == resname(i)) SeqID = i
    end do

  end function SeqID

  !below function calculates SHD value for a given sequence
  function getSHD(s) result(y)

    implicit none

    character(len=*), intent(in) :: s
    real(8) :: y
    integer :: i, j, n
    real(8) :: x, xt, wt

    y = 0.0
    n = len_trim(s)
    do i = 2, n
       do j = 1, i-1
          x = dble(i-j)
          xt = x**(-0.5)
          wt = 1.0 ! exp(-x*x/25.0)
          y = y + hpres(SeqID(s(i:i)),SeqID(s(j:j)))*xt*wt
       end do
    end do 
    y = y/dble(n)
   
  end function getSHD

  !below function calculates SCD value for a given sequence
  function getSCD(s) result(y)

    implicit none

    character(len=*), intent(in) :: s
    real(8) :: y
    integer :: i, j, n
    real(8) :: x, xt, wt

    y = 0.0
    n = len_trim(s)
    do i = 2, n
       do j = 1, i-1
          x = dble(i-j)
          xt = x**(0.5)
          wt = 1.0 ! exp(-x*x/25.0)
          y = y + qres(SeqID(s(i:i)))*qres(SeqID(s(j:j)))*xt*wt
       end do
    end do 
    y = y/dble(n)
   
  end function getSCD

  !below function computes the number of arguments in a given line of a file 
  function LineArgs(line) result(n)

    implicit none

    character(len=*), intent(in) :: line
    integer :: n
    character(len=len(line)) :: txt
    integer :: i1

    n = 0
    i1 = 1
    call read_lineblock(i1, line, txt)
    do while (i1 <= len(line) .and. len_trim(txt) > 0 .and. &
       txt(1:1) /= '#')
       n = n + 1
       call read_lineblock(i1, line, txt)
    end do

  end function LineArgs

  !below subroutines gets the actual arguments in a given line of a file
  subroutine GetLineArgs(line, n, df)

    implicit none

    character(len=*), intent(in) :: line
    integer, intent(in) :: n
    character(len=255), intent(out) :: df(n)
    character(len=len(line)) :: txt
    integer :: i1, nt

    df = ''
    nt = 0
    i1 = 1
    call read_lineblock(i1, line, txt)
    do while (i1 <= len(line) .and. len_trim(txt) > 0 .and. &
          txt(1:1) /= '#')
       nt = nt + 1
       if (nt > n) stop 'Number of fields exceed the limit!'
       df(nt) = trim(txt)
       call read_lineblock(i1, line, txt)
    end do

  end subroutine GetLineArgs

  subroutine read_lineblock(i, line, bl)

    implicit none

    integer, intent(inout) :: i
    character(len=*), intent(in) :: line
    character(len=len(line)), intent(out) :: bl
    integer :: i1, i2

    bl = ''

    i1 = i
    do while (line(i1:i1) == ' ' .and. i1 < len(line))
       i1 = i1 + 1
    end do
    i2 = i1
    do while (line(i2:i2) /= ' ' .and. i2 < len(line))
       i2 = i2 + 1
    end do
    read (line(i1:i2), '(a)') bl
    i = min(i2, len(line))

  end subroutine read_lineblock 

end program main
