program main

  integer                       :: i, nt, nsteps, dumping_freq
  character(200)                :: initialFname, cmd_name, fname, hiddenInput
  character(200)                :: inputFolder, controlFile, hiddenOutput
  character(200)                :: snap_name, outputFolder, controlOutput
  character(200)                :: pname, sensorFile, infoFile
  character(200)                :: fullstate_controlFile
  character(20)                 :: buffer
  character(8)                  :: stepString
  logical                       :: here

  ! grid variables
  real(8), dimension(249,249)   :: x, y
  real(8), dimension(62001)     :: actuator, Kgain
  real(8)                       :: xc1, xc2, yc1, yc2, a, c, r1, r2
  integer                       :: nx, ny, j
  real(8)                       :: total
  real(8), dimension(:), allocatable  :: omega1

  ! random seed variables
  integer, dimension(:), allocatable :: seed
  integer                            :: clock, nseed, inputseed

  ! ROM variables
  integer                              :: nhidden, nobserve
  real(8), dimension(:,:), allocatable :: Am, B, Cm, L, K
  real(8), dimension(:), allocatable   :: xhidden, yobserve
  real(8)                              :: strength, full_strength

  ! sensor variables
  real(8)                              :: xs1, xs2, ys1, ys2, sigma
  real(8), dimension(:,:), allocatable :: sensorC
  integer                              :: iter, nsensors

  print*, 'enter simulation name:'
  read*, pname

  print*, 'enter the number of iterations for actuator: '
  read*, nt

  print*, 'timestep of controller: '
  read*, nsteps

  print*, 'starting iteration: '
  read*, istart

  print*, 'dumping frequency: '
  read*, dumping_freq

  print*, 'seed for the random number generator: '
  read*, inputseed

  print*, 'enter input files folder name: '
  read*, inputFolder

  print*, 'enter the sensor filename: '
  read*, sensorFile

  print*, 'enter the controller file name: '
  read*, controlFile

  print*, 'enter fullstate control file name:'
  read*, fullstate_controlFile

  print*, 'enter the initial condition restart file name: '
  read*, initialFname

  print*, 'enter the folder name for snapshots: '
  read*, outputFolder

  print*, 'enter the control output file name: '
  read*, controlOutput

  print*, 'enter hidden state input file name: '
  read*, hiddenInput

  print*, 'enter hidden state output prefix: '
  read*, hiddenOutput

  call execute_command_line('mkdir -p '//adjustl(trim(outputFolder)))

  infoFile = trim(outputFolder)//'/'//trim(adjustl(pname))//'_info.txt'
  inquire(file=infoFile, exist=here)
  if(here) then
    print*, 'info file already exists, ', infoFile
    call exit(1)
  else
    open(101,file=infoFile,action='write',status='new',form='formatted')
  endif

  write(101,*) 'project name --- ', pname
  write(101,*) 'iterations --- ', nt
  write(101,*) 'step size --- ', nsteps
  write(101,*) 'start index --- ', istart
  write(101,*) 'dumping frequency --- ', dumping_freq
  write(101,*) 'input seed --- ', inputseed
  write(101,*) 'input folder --- ', inputFolder
  write(101,*) 'sensor file --- ', sensorFile
  write(101,*) 'control file --- ', controlFile
  write(101,*) 'fullstate control file --- ', fullstate_controlFile
  write(101,*) 'input restart file --- ', initialFname
  write(101,*) 'output folder --- ', outputFolder
  write(101,*) 'control output --- ', controlOutput
  write(101,*) 'hidden input --- ', hiddenInput
  write(101,*) 'hidden output --- ', hiddenOutput

  close(101)

  ! random seed
  call random_seed(size = nseed)
  allocate(seed(nseed))
  call system_clock(count = clock)
  if (inputseed < 0) then
    seed = clock + 37 * (/ (i - 1, i = 1, nseed) /)
  else
    seed = inputseed
  endif

  call random_seed(put = seed)

  cmd_name = trim('ln -s ' // trim(inputFolder) // '/mlDMDcontrol_01.cholesky ' &
                  // trim(pname) // '_01.cholesky')
  call execute_command_line(cmd_name)

  cmd_name = trim('ln -s ' // trim(inputFolder) // '/mlDMDcontrol_02.cholesky ' &
                  // trim(pname) // '_02.cholesky')
  call execute_command_line(cmd_name)

  cmd_name = trim('ln -s ' // trim(inputFolder) // '/mlDMDcontrol_03.cholesky ' &
                  // trim(pname) // '_03.cholesky')
  call execute_command_line(cmd_name)

  initialFname = trim(inputFolder)//'/'//trim(initialFname)
  controlOutput = trim(outputFolder)//'/'//trim(controlOutput)
  sensorFile = trim(inputFolder)//'/'//trim(sensorFile)
  controlfile = trim(inputfolder)//'/'//trim(controlfile)
  fullstate_controlFile = trim(inputfolder)//'/'//trim(fullstate_controlFile)
  hiddenInput = trim(inputFolder)//'/'//trim(hiddenInput)
  hiddenOutput = trim(outputFolder)//'/'//trim(hiddenOutput)

  print*, 'project name is: ', pname
  print*, 'initial restart files is: ', initialFname
  print*, 'the actuator output file is: ', controlOutput
  print*, 'control input file name is: ',controlFile
  print*, 'hidden variable input file name is: ',hiddenInput
  print*, 'hidden variable output file name is: ',hiddenOutput

  ! cp inputfile into <pname>_ic.bin

  allocate(omega1(62001))
  call read_file(initialFname, omega1)
  write(stepString, '(I5.5)') nsteps
  cmd_name = trim('cp '//trim(initialFname) &
                  //' '//trim(pname)//trim(adjustl(stepString))//'.bin')
  print*, trim(cmd_name)
  call execute_command_line(cmd_name)
  fname = trim(pname)//'_ic.bin'
  call write_file(nsteps, pname, omega1, fname)

  deallocate(omega1)

  ! read x and y values
  open(101, file=trim(inputFolder)//'/grid.dat', form='formatted')

  do j=1,249
    do i=1,249
      read(101,*) x(i,j)
    enddo
  enddo

  do j=1,249
    do i=1,249
      read(101,*) y(i,j)
    enddo
  enddo

  close(101)

  print*, maxval(x), minval(x), maxval(y), minval(y)

  ! impulse parameters
  c = 2.0d0
  a = 20.0d0
  xc1 = 0.0d0
  xc2 = 0.0d0
  yc1 = 1.3423d0
  yc2 = 0.89d0

  actuator = 0.d0

  nx = 250
  ny = 250

  do j=1,ny-1
    do i=1,nx-1
      r1 = sqrt((x(i,j) - xc1)**2 + (y(i,j) - yc1)**2)
      r2 = sqrt((x(i,j) - xc2)**2 + (y(i,j) - yc2)**2)
      actuator(i + (nx-1)*(j-1)) = actuator(i + (nx-1)*(j-1)) + &
                                  c*(1 - a*r1**2)*exp(-a*r1**2) - &
                                  c*(1 - a*r2**2)*exp(-a*r2**2)
    enddo
  enddo

  print*, maxval(actuator), minval(actuator)

  ! sensor parameters
  sigma = 0.1d0
  nx = 250
  ny = 250

  open(101, file = trim(adjustl(sensorFile)), form = 'formatted')

  read(101, *) nsensors

  allocate(sensorC(nsensors, 62001))
  sensorC = 0.d0

  do iter = 1, nsensors

    read(101, *) xs1
    read(101, *) ys1


	  do j=1,ny-1
	    do i=1,nx-1
	      r1 = sqrt((x(i,j) - xs1)**2 + (y(i,j) - ys1)**2)
!	      r2 = sqrt((x(i,j) - xs2)**2 + (y(i,j) - ys2)**2)
	      sensorC(iter, i + (nx-1)*(j-1)) = sensorC(iter, i + (nx-1)*(j-1)) + &
	                                     exp(-r1**2/2/sigma**2)
!	      sensorC(2, i + (nx-1)*(j-1)) = sensorC(2, i + (nx-1)*(j-1)) + &
!	                                     exp(-r2**2/2/sigma**2)
	    enddo
	  enddo

  enddo

  close(101)


  ! read the file containig the feedback gain matrix
!  open(101, file=trim(adjustl(fullstate_controlFile)), form='formatted')
!
!  do i=1,62001
!    read(101,*) Kgain(i)
!  enddo
!
!  close(101)


  ! initialize the reduced order model
  open(101, file = trim(adjustl(controlFile)), form='formatted')

  read(101, *) nhidden
  read(101, *) nobserve

  allocate(Am(nhidden,nhidden), B(nhidden,1), Cm(nobserve,nhidden))
  allocate(L(nhidden,nobserve), K(1,nhidden))
  allocate(xhidden(nhidden), yobserve(nobserve))

  do i = 1,nhidden
    do j = 1,nhidden
      read(101, *) Am(j,i)
    enddo
  enddo

  do i = 1,nhidden
    read(101, *) B(i,1)
  enddo

  do i = 1,nhidden
    do j = 1,nobserve
      read(101, *) Cm(j,i)
     enddo
   enddo

   do i = 1,nobserve
    do j = 1,nhidden
      read(101, *) L(j,i)
    enddo
  enddo

  do i = 1,nhidden
    read(101, *) K(1,i)
  enddo

  close(101)

  if ( nsensors.ne.nobserve ) then
    print*, 'number of sensors do not match'
    call exit(1)
  endif

  call execute_command_line('rm -rf '//trim(outputFolder)//'/mlDMDcontrol.force')
  call execute_command_line('rm -rf '//trim(outputFolder)//'/mlDMDcontrol.energy')

  strength = 0.d0
  total = 0.d0

  if (istart == 0) then
    xhidden = 0.d0
    yobserve = 0.d0
  else
    call readHiddenInput(nhidden, xhidden, hiddenInput)
  endif


  print*, 'number of sensors is ',nsensors
  print*, 'number of observers is ', nobserve

  ! loop for controlled simulation
  do i=istart+1,istart+nt

    ! advance the simulation simulation
    call advanceOmega(nsteps, pname)

    call execute_command_line('cat '//trim(pname)//'.force >> ' &
                              //trim(adjustl(outputFolder))//'/forceFile.txt')
    call execute_command_line('cat '//trim(pname)//'.energy >> ' &
                              //trim(adjustl(outputFolder))//'/energyFile.txt')

    ! get observation
    call getObservation(nsteps, pname, nobserve, yobserve, total, &
                        sensorC, inputFolder)

    ! get new actuator value by advancing estimate
    call getControl(nhidden, nobserve, &
                    xhidden, yobserve, Am, B, Cm, L, K, strength)

    !call getFullStrength(pname, Kgain, full_strength, inputFolder)

    ! write the strength value
    inquire(file=trim(adjustl(controlOutput)), exist=here)
    if(here) then
      open(100,file=trim(adjustl(controlOutput)),position='append', &
                 action='write', status='old')
    else
      open(100,file=trim(adjustl(controlOutput)),action='write',status='new')
    endif

    write(100,*) i, strength, total, yobserve(1)!, full_strength

    close(100)

    print*, 'the applied strength and observation is ', strength, yobserve(1), rnorm()

    if (isnan(yobserve(1))) stop 'solution is nan'


    ! apply controller
    call applyControl(nsteps, pname, actuator, strength, inputFolder)
    !call applyControlObservation(nsteps, nobserve, yobserve, &
    !                             actuator, Cm, strength, inputFolder)

    write(buffer,'(i5.5)') i
    snap_name = trim(adjustl(trim(outputFolder)//'/snapShot_'//trim(pname)//trim(buffer)//'.plt'))
    cmd_name = trim('cp '//trim(pname)//'00000.plt '//trim(snap_name))
    if ( mod(i,dumping_freq) == 0 ) then
      call execute_command_line(cmd_name)
      call execute_command_line('cp '//trim(pname)//'_ic.bin '// &
                        trim(outputFolder)//'/restart_'//trim(pname)//trim(buffer)//'.bin')
      call writeHiddenOutput(i, nhidden, xhidden, hiddenOutput)
    endif



  enddo

  ! write everything before terminating the code
  write(buffer,'(i5.5)') i-1
  snap_name = trim(adjustl(trim(outputFolder)//'/snapShot_'//trim(pname)//trim(buffer)//'.plt'))
  cmd_name = trim('cp '//trim(pname)//'00000.plt '//trim(snap_name))
  if ( mod(i,dumping_freq) .ne. 0 ) then
    call execute_command_line(cmd_name)
    call execute_command_line('cp '//trim(pname)//'_ic.bin '// &
                    trim(outputFolder)//'/restart_'//trim(pname)//trim(buffer)//'.bin')
    call writeHiddenOutput(i-1, nhidden, xhidden, hiddenOutput)
  endif

  deallocate(Am,B,Cm,L,K,xhidden,yobserve)

end program main

subroutine getControl(nhidden, nobserve, xhidden, yobserve, A, B, C, L, K, strength)

  integer                               :: nhidden, nobserve
  real(8), dimension(nhidden, nhidden)  :: A
  real(8), dimension(nhidden, 1)        :: B
  real(8), dimension(nobserve, nhidden) :: C
  real(8), dimension(nhidden, nobserve) :: L
  real(8), dimension(1, nhidden)        :: K

  ! x = Ax + Bu + L(y - C(Ax + Bu))
  ! u = -Kx
  real(8), dimension(nhidden)        :: xhidden, xtemp
  real(8), dimension(nobserve)       :: yobserve, ytemp
  real(8)                            :: strength

  integer                            :: i,j

  xtemp = 0.d0

  ! xtemp = Ax
  do i = 1,nhidden
    do j = 1,nhidden
      xtemp(i) = xtemp(i) + A(i,j) * xhidden(j)
    enddo
  enddo

  ! xtemp = Ax + Bu
  xtemp(1:nhidden) = xtemp(1:nhidden) + B(1:nhidden,1) * strength

  ytemp = 0.d0

  ! ytemp = C(Ax + Bu)
  do i = 1,nobserve
    do j = 1,nhidden
      ytemp(i) = ytemp(i) + C(i,j) * xtemp(j)
    enddo
  enddo

  ! xtemp = Ax + Bu + L(y - C(Ax + Bu))
  do i = 1,nhidden
    do j = 1,nobserve
      xtemp(i) = xtemp(i) + L(i,j) * (yobserve(j) - ytemp(j))
    enddo
  enddo

  strength = 0.d0

  ! strength = -Kxtemp
  do i = 1,nhidden
    strength = strength - K(1,i) * xtemp(i)
  enddo

  ! x = xtemp
  xhidden = xtemp

end subroutine getControl

subroutine advanceOmega(nsteps, pname)

  integer                         :: nsteps
  character(200)                  :: pname
  character(400)                  :: cmdname
  character(8)                    :: stepString

  write(stepString, '(I5)') nsteps

  cmdname = './ibpm -name '//trim(adjustl(pname))//' -tecplot '//trim(adjustl(stepString)) &
            //' -restart '//trim(adjustl(stepString)) &
            //' -nx 250 -ny 250 -ngrid 5 -length 5 -xoffset -2 -yoffset -2.5 ' &
            //'-energy 1 -geom flatplate35.geom -model nonlinear' &
            //' -ic '//trim(pname)//'_ic.bin -dt 0.01 -nsteps ' &
            //trim(adjustl(stepString))//' > '//trim(pname)//'_log.txt'

  call execute_command_line(cmdname)
  print*, cmdname

end subroutine advanceOmega

subroutine getObservation(nsteps, pname, nobserve, yobserve, total, &
                          C, inputFolder)

  integer                           :: nobserve, nsteps, i
  real(8), dimension(nobserve)      :: yobserve
  character(8)                      :: stepString
  character(200)                    :: pname, inputFolder
  character(200)                    :: fileName, fileMean, fname
  real(8), dimension(62001)         :: yFull, yOut, yMean, actuator
  real(8), dimension(nobserve,62001):: C
  real(8)                           :: strength, total

  write(stepString, '(I5.5)') nsteps

  fileName = trim(pname)//trim(adjustl(stepString))//'.bin'

  print*, 'reading from ', fileName

  fileMean = trim(adjustl(inputFolder))//'/steadystate.bin'

  call read_file(fileName, yFull)

  call read_file(fileMean, yMean)

  ! apply noise and make observation
  total = 0.d0
  do i = 1,62001
    ! system noise
    yOut(i) = yFull(i) + 0.01d0 * rnorm()
    total = total + (yOut(i) - yMean(i))**2
  enddo

  do i = 1, nobserve
    yobserve(i) = 0.d0
    do j = 1,62001
      yobserve(i) = yobserve(i) + C(i,j) * (yOut(j) - yMean(j))
    enddo
    ! observation error
    yobserve(i) = yobserve(i) + 0.01d0 * rnorm()
  enddo

  ! write the restart file
  fname = trim(pname)//'_temp.bin'
  call write_file(nsteps, pname, yOut, fname)

end subroutine getObservation


subroutine applyControl(nsteps, pname, actuator, strength, inputFolder)

  integer                           :: nsteps, i
  character(200)                    :: pname, inputFolder
  character(200)                    :: fileName, fname
  real(8), dimension(62001)         :: yFull, yOut, actuator
  real(8)                           :: strength

  fileName = trim(pname)//'_temp.bin'

  call read_file(fileName, yFull)

  ! apply actuation
  yOut = yFull + actuator * strength


  ! write the restart file
  fname = trim(pname)//'_ic.bin'
  call write_file(nsteps, pname, yOut, fname)

end subroutine applyControl

subroutine getFullStrength(pname, Kgain, full_strength, inputFolder)

  integer                           :: i
  character(200)                    :: pname, inputFolder
  character(200)                    :: fileName
  real(8), dimension(62001)         :: Kgain, yFull, yMean
  real(8)                           :: full_strength

  fileName = trim(pname)//'_temp.bin'
  call read_file(fileName, yFull)

  fileName = trim(adjustl(inputFolder))//'/steadystate.bin'
  call read_file(fileName, yMean)

  full_strength = 0.d0
  do i = 1,62001
    full_strength = full_strength - Kgain(i)*(yFull(i) - yMean(i))
  enddo

end subroutine getFullStrength

subroutine readHiddenInput(nhidden, xhidden, hiddenInput)

  integer                         :: i, nhidden
  real(8), dimension(nhidden)     :: xhidden
  character(200)                  :: hiddenInput

  open(101, file=trim(adjustl(hiddenInput)), form = 'formatted')

  do i = 1,nhidden
    read(101,*) xhidden(i)
  enddo

  close(101)

end subroutine readHiddenInput

subroutine writeHiddenOutput(iter, nhidden, xhidden, hiddenOutput)

  integer                           :: iter, i, nhidden
  real(8), dimension(nhidden)       :: xhidden
  character(200)                    :: hiddenOutput
  character(20)                     :: buffer

  write(buffer,'(i5.5)') iter
  open(101, file=trim(adjustl(hiddenOutput))//'_'//&
                 trim(adjustl(buffer))//'.txt', action='write')

  do i = 1,nhidden
    write(101,*) xhidden(i)
  enddo

  close(101)

end subroutine writeHiddenOutput
