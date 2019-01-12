! smc_driver_mpi.f90 -- MPI driver for Sequential Monte Carlo Algorithm
!   Ed Herbst [edward.p.herbst@frb.gov]

program smc_driver
  use f95_precision, only: wp => dp

  use mcmc

  !-----------------------------------------------------------------------------
  ! model selection
  ! use {nkmp,sw,rbc,bimodal_dsge,news,edo,edo_exp}
  !-----------------------------------------------------------------------------
  use sw !bimodal_dsge
  use mpi

  implicit none

  !include 'mpif.h'

  integer :: mpierror, rank, nproc
  integer :: npart, nintmh, nphi, irep, nT, nblocks, nomp, ngap
  logical :: init_prior, do_geweke, phi_bend

  integer, allocatable :: Tvec(:)
  real(wp), allocatable :: phi(:), wtsim(:), parasim(:,:)
  real(wp) :: kapp, bcoeff
  real(wp) :: lik0, p2(2, npara)

  character(len=144) :: basedir
  character(len=144) :: cnpart, cnintmh, cnphi, ckap, cnblocks
  character(len=144) :: cinitdist, cirep, cmethod
  character(len=144) :: initpart, initw
  character(len=144) :: arg, fstr

  integer :: i, j, ii, m, mm !RECA
  integer, allocatable :: break_points(:)
  character(len=244) :: chareca ! RECA


  ! mask stuff
  integer :: neffpara
  integer, allocatable :: pmaskinv(:), ind(:)
  real(wp), allocatable :: eye(:,:)

  ! file i/o
  logical :: direxists
  character(len=244) :: phifile, parafile, wghtfile, llikfile, postfile, zifile, rfile, essfile

  character(len=144) :: hyp_str, hyp_fstr
  logical :: fixed_hyper

  ! random numbers
  type (VSL_STREAM_STATE) :: stream
  integer :: brng, seed, methodu, method, mstorage, time_array(8), errcode

  real(wp), allocatable :: uu(:), u(:), eps(:,:)
  real(wp) :: mu(npara), var(npara, npara), scale, ahat, zero_vec(npara)
  integer, allocatable :: paraind(:)

  ! containers
  real(wp), allocatable :: arate(:,:), loglh(:), prio(:), ESS(:), zi(:)
  real(wp), allocatable :: incwt(:), wtsq(:)
  integer, allocatable :: acptsim(:,:), resamp(:)

  ! RECA
  real(wp), allocatable :: save_mixr(:,:), save_step_lik(:,:), save_step_pr(:,:)
  real(wp), allocatable :: nodemixr(:,:), nodesteplik(:,:), nodesteppr(:,:)
  !!!!!!!

  real(wp), allocatable :: nodepara(:,:), nodeloglh(:), nodewt(:), nodeprio(:)
  real(wp), allocatable :: nodeeps(:,:), nodeu(:), nodevar(:,:)
  integer, allocatable :: nodeacpt(:,:)

  logical :: do_resample

  !------------------------------------------------------------
  ! Initialize MPI
  !------------------------------------------------------------
  call mpi_init(mpierror)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)


  ! base directory for results
  call get_environment_variable('OUTPUT_BASE', basedir)
  basedir = trim(adjustl(basedir))

  !------------------------------------------------------------
  ! default options
  !------------------------------------------------------------
  npart   = 12000               ! number of particles

  ! tempering schedule
  cmethod = 'smc'               ! character holder
  do_geweke = .false.           ! FALSE = simulated tempering SMC
  nphi    = 500                 ! number of heating stages
  phi_bend = .true.             ! do bending
  bcoeff  = 2.1_wp              ! lambda coefficient for heating

  ! MCMC parameters
  nintmh  = 1                   ! number of intermediate MH steps per stage
  nblocks = 1                   ! number of blocks !667
  kapp    = 0.9_wp              ! weight on standard RWMH proposal
  ckap    = 'mix-'              ! character holder
  scale = 0.4_wp                ! starting scale on variance


  ! initialization
  init_prior = .true.           ! initialize particles by sampling prior
  cinitdist = 'prior'           ! character holder

  ! hyperparameters
  hyp_str = ''
  hyp_fstr = ''
  fixed_hyper = .false.

  ! only one of these is used depending on SMC type
  nT = nobs
  initpart = initfile           ! initial particles (if init_prior = .false.)
  initw = initwt                ! initial weight (if init_prior = .false.)

  call omp_set_num_threads(12)  ! OMP threads
  call mkl_set_num_threads(1)   ! MKL threads


  irep = 1                      ! run index

  !------------------------------------------------------------
  ! read in command line options
  !------------------------------------------------------------
  do i = 1, command_argument_count()
     call get_command_argument(i, arg)

     select case(arg)
     case('-p', '--nphi')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') nphi
     case('-n', '--npart')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') npart
     case('-m', '--nintmh')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') nintmh
     case('-t', '--studt')
        init_prior = .false.
        cinitdist = 'studt'
     case('-pr', '--prior')
        init_prior = .true.
        cinitdist = 'prior'
     case('-b', '--bend')
        phi_bend = .true.
        call get_command_argument(i+1, arg)
        read(arg, '(f)') bcoeff
     case('-u','--proposal rwmh')
        kapp = 1.0_wp
        ckap = ''
     case('-o', '--nblocks')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') nblocks
     case('-g', '--geweke')
        do_geweke = .true.
        cmethod = 'geweke'
     case('-i', '--initrep')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') irep
     case('--nomp')
        call get_command_argument(i+1, arg)
        read(arg, '(i)') nomp
        call omp_set_num_threads(nomp)
     case('--diffuse-kalman')
        ! not yet implemtn

     case('--output-dir')
        call get_command_argument(i+1, arg)
        basedir = arg
     case('-h','--help')
        call print_help()
        stop
     case('--fixed-hyper')
        call get_command_argument(i+1, hyp_fstr)
        hyp_str = '-hypfixed'
        fixed_hyper = .true.
     case('-f', '--init-files')
!        call get_command_argument(i+1, initpart)
!        call get_command_argument(i+2, initwt)
     end select
  end do

  !-------------------------------------------------------------------------
  ! Create output directory
  !-------------------------------------------------------------------------
  write(cnphi, '(i4)') nphi
  write(cnintmh, '(i2)') nintmh
  write(cnpart, '(i8)') npart
  write(cirep, '(i2)') irep
  write(cnblocks, '(i2)') nblocks

  fstr = trim(adjustl(cmethod))//'-'//trim(mname)//'-'//trim(ckap)//'npart-' &
       //trim(adjustl(cnpart))//'-nintmh-'//trim(adjustl(cnintmh))//'-nphi-' &
       //trim(adjustl(cnphi))//'-'//trim(cinitdist)//'-b'//trim(adjustl(cnblocks))&
       //trim(adjustl(hyp_str))//'-trial'//trim(adjustl(cirep))

  if (phi_bend == .true.) then
     fstr = trim(adjustl(fstr))//'-phibend'
  endif

  if (rank == 0) then
     print*,'Writing output to...'
     write(*,'(A)'),fstr

     print*,'# processors = ', nproc

     if (not(mod(npart, nproc) == 0)) then
        print*,'ERROR npart must be divisible by nsim!'
        stop
     end if
  end if

  !-------------------------------------------------------------------------
  ! read in data
  !-------------------------------------------------------------------------
  call read_in_from_files()

  !-------------------------------------------------------------------------
  ! create heating schedule and check mod(npart/nproc) == 0
  !  see the middle of page 9
  !-------------------------------------------------------------------------
  allocate(phi(nphi))
  do i = 1, nphi
     phi(i) = (i-1.0_wp)/(nphi - 1.0_wp)
  end do

  if (phi_bend == .true.) then
     if (rank == 0) then
        print*,'bending phi'
     end if
     phi = phi**bcoeff
  end if

  !-------------------------------------------------------------------------
  ! allocate memory on nodes
  !-------------------------------------------------------------------------
  ngap = int(npart*1.0_wp/nproc)
  allocate(nodepara(npara,ngap), nodeloglh(ngap), nodewt(ngap), &
       nodeprio(ngap), nodeacpt(npara,ngap), nodevar(npara,npara))
  nodeacpt = 0
  allocate(ESS(nphi), zi(nphi))

  ! RECA - MIX
  allocate(nodemixr(nblocks*nintmh, ngap))
  allocate(nodesteplik(nblocks*nintmh, ngap))
  allocate(nodesteppr(nblocks*nintmh*npara, ngap))
  !-------------------------------------------------------------------------
  ! draw initial particles
  !-------------------------------------------------------------------------
  if (rank == 0) then
     allocate(parasim(npara, npart), wtsim(npart), loglh(npart),&
     & prio(npart), incwt(npart), wtsq(npart))
     allocate(paraind(npart), acptsim(npara,npart),resamp(nphi))

     allocate(save_mixr(nblocks*nintmh, npart))
     allocate(save_step_lik(nblocks*nintmh, npart))
     allocate(save_step_pr(nblocks*nintmh*npara, npart))

     parasim = transpose(priorrand(npart, pshape, pmean, pstdd, pmask, pfix))

     write(*,'(a$)') 'Evaluating initial likelihoods. '

  endif


  call mpi_barrier(MPI_COMM_WORLD, mpierror)
  call scatter(everything=.false.)
  call mpi_barrier(MPI_COMM_WORLD, mpierror)

  call initial_draw_from_prior(nodeloglh, nodeprio, nodepara, nodewt, ngap, npara)
  call mpi_barrier(MPI_COMM_WORLD, mpierror)
  call gather(everything=.false.)
  call mpi_barrier(MPI_COMM_WORLD, mpierror)


  fstr = trim(basedir)//trim(fstr)
  !hyp_fstr = trim(basedir)//trim(hyp_fstr)
  if (rank == 0) then
     print*, 'Saving files to:', fstr
     inquire(directory=fstr, exist=direxists)
     if (direxists) then
        print *, 'directory already exists'
     else
        call system('mkdir '//fstr)
     endif
     phifile = trim(fstr)//'/phi.txt'

     open(1, file=phifile, action='write')
     do i = 1, nphi
        write(1, '(1f)') phi(i)
     end do
     close(1)
  endif



  !------------------------------------------------------------
  ! construct mask arrary (refactor this)
  !------------------------------------------------------------
  neffpara = npara - sum(pmask)

  print*, neffpara
  print*, npara

  allocate(pmaskinv(neffpara), ind(neffpara), eye(npara, npara), eps(npart*nintmh,npara))
  allocate(nodeeps(ngap*nintmh,npara), break_points(nblocks+1))

  j = 1
  do i = 1, npara
     eye(i,i) = 1.0_wp
     if (pmask(i) == 0) then
        pmaskinv(j) = i
        j = j + 1
     end if
  end do

  print *, pmaskinv

  call date_and_time(values=time_array)


  !----------------------------------------------------------------------------
  ! Initialize random deviates only using master
  !---------------------------------------------------------------------------
  ! for generating random numbers
  eye = 0.0_wp
  zero_vec = 0.0_wp

  if (rank == 0) then
     brng = VSL_BRNG_MT19937
     seed = 2105 !mod(sum(time_array),10000)
     print*,'RANDOM SEED=',seed
     methodu=VSL_RNG_METHOD_UNIFORM_STD
     method=VSL_RNG_METHOD_GAUSSIANMV_ICDF
     mstorage=VSL_MATRIX_STORAGE_FULL

     allocate(uu(nphi*npart), u(npart*nphi*nintmh*nblocks))
     errcode = vslnewstream(stream, brng, seed)
     errcode = vdrnguniform(methodu, stream, nphi*npart, uu, 0.0_wp, 1.0_wp)
     errcode = vdrnguniform(methodu, stream, npart*nphi*nintmh*nblocks, u, 0.0_wp, 1.0_wp)

     allocate(arate(nphi, npara))
     arate(:,:) = 0.0_wp

     call write_model_para(fstr)

     ! RECA PRINTING
     open (unit=69, file = trim(fstr)//'/'//'stepprobs.txt', action='write')
     do m = 1, npart*nphi*nintmh*nblocks
        write (69, '(f20.12)') u(m)
     end do
     close(69)
  else
     allocate(uu(nphi*npart), u(npart*nphi*nintmh*nblocks))
  endif


  allocate(nodeu(ngap*nintmh*nphi*nblocks))
  call mpi_scatter(u, ngap*nphi*nintmh*nblocks, MPI_DOUBLE_PRECISION, nodeu, ngap*nphi*nintmh*nblocks, &
       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)

  !-----------------------------------------------------------------------------
  ! INITIALIZE WEIGHTS
  !-----------------------------------------------------------------------------
  if (rank == 0) then
     wtsim = 1.0_wp/npart
     wtsim = wtsim / sum(wtsim)

     ESS(1) = 1.0_wp / sum(wtsim**2)
     zi(1) = 1.0_wp

     call write_files(1)

     zifile = trim(fstr)//'/zratios.txt'
     open(3, file=zifile, action='write')
     write(3, '(f)') zi(1)

     resamp = 0

     ! load the resampling scheme if true
     if (fixed_hyper==.true.) then
        rfile = trim(hyp_fstr)//'/'//'resamp.txt'
        open(1,file=rfile,action='read')
        do i = 1,nphi
           read(1,*) resamp(i)
        end do
        close(1)
     end if
  end if

  !-----------------------------------------------------------------------------
  ! BEGIN SEQUENTIAL MONTE CARLO
  !-----------------------------------------------------------------------------
  do i = 2, nphi

     if (rank == 0) then

        !---------------------------------------------------------------------------
        ! [C]orrection
        ! This is in importance sampling step
        !---------------------------------------------------------------------------

        !------------------------------------------------------------
        ! incremental weight = exp((phi[i]-phi[i-1])*loglh)
        !------------------------------------------------------------
        call vdexp(npart, (phi(i) - phi(i-1))*loglh, incwt)

        ! total weight
        wtsim = wtsim * incwt

        ! incremental MDD
        zi(i) = sum(wtsim)
        write(3, '(f)') zi(i)

        ! normalize weights
        wtsim = wtsim / zi(i)


        !----------------------------------------------------------------------------
        ! [S]election
        ! (Potentially) Resample the particles
        !----------------------------------------------------------------------------

        ! compute effective sample size
        call vdsqr(npart, wtsim, wtsq)
        ESS(i) = 1.0_wp/sum(wtsq)

        write(*,'(a,f4.2,a,f$)') 'phi[', phi(i), '] -- [ESS = ', ESS(i)


        if (fixed_hyper==.true.) then
           do_resample = resamp(i)==1
        else
           do_resample = ESS(i) < npart/2.0_wp
        end if

        if (do_resample) then
           ! do resampling step, if necessary
           write (*, '(a$)') ' resample]'
           call multinomial_resampling(npart, wtsim, uu((i-1)*npart+1:i*npart), paraind)

           ! RECA PRINTING
           write (chareca,'(i3.3)') i
           open (unit=69, file = trim(fstr)//'/'//trim(chareca)//'resamp.txt', action='write')
           do m = 1, npart
              write (69, '(100i)') paraind(m)
           end do
           close(69)

           ! resample
           parasim = parasim(:, paraind)
           loglh = loglh(paraind)
           prio = prio(paraind)

           wtsim = 1.0_wp / npart

           resamp(i) = 1

        else
           write (*, '(a$)') '         ]'
        end if


        !-------------------------------------------------------------------------------
        ! [M]utation
        ! Propagate the particles to towards areas of high density
        !  (potentially using adaptive hyperparameters).
        !
        ! Note that this step continues behind the hanging 'end if' which exists
        !  because of MPI.
        !-------------------------------------------------------------------------------

        ! find mu and variance for MH step
        call compute_mean_and_variance(mu, var, parasim, wtsim, npara, npart)

        ! generate random blocks
        ind = pmaskinv
        call generate_random_blocks(0.0_wp, neffpara, nblocks, ind, break_points)

        if (fixed_hyper==.true.) then
           call load_hyperparameters(i)
        end if

        ! save output
        call write_hyperparameters(i)

        errcode=vdrnggaussian(method,stream,npart*nintmh*npara,eps,0.0_wp,scale)
        acptsim = 0


     endif

     ! scatter info to nodes
     call scatter(everything=.true.)
     call mpi_barrier(MPI_COMM_WORLD, mpierror)

     ! do particle mutation
     call transition_kernel(nodeloglh, nodeprio, nodepara, nodeacpt, nodeeps, ngap, var, nodemixr, nodesteppr, nodesteplik)

     ! send info back to master node
     call mpi_barrier(MPI_COMM_WORLD, mpierror)
     call gather(everything=.true.)
     call mpi_barrier(MPI_COMM_WORLD, mpierror)

     ! RECA PRINTING mixr DEETS
     if (rank == 0) then
        write (chareca,'(i3.3)') i
        open (unit=69, file = trim(fstr)//'/'//trim(chareca)//'mixr.txt', action='write')
        do m = 1, npart
           do mm = 1, nblocks*nintmh
              if (mm == nblocks*nintmh) then
                 write (69, '(100f)') save_mixr(mm, m)
              else
                 write (69, '(100f)', advance='no') save_mixr(mm, m)
              end if
           end do
        end do
        close (69)
        open (unit=69, file = trim(fstr)//'/'//trim(chareca)//'step_lik.txt', action='write')
        do m = 1, npart
           do mm = 1, nblocks*nintmh
              if (mm == nblocks*nintmh) then
                 write (69, '(100f)') save_step_lik(mm, m)
              else
                 write (69, '(100f)', advance='no') save_step_lik(mm, m)
              end if
           end do
        end do
        close (69)
        open (unit=69, file = trim(fstr)//'/'//trim(chareca)//'step_pr.txt', action='write')
        do m = 1, npart
           do mm = 1, nblocks*nintmh*npara
              if (mm == nblocks*nintmh*npara) then
                 write (69, '(100f)') save_step_pr(mm, m)
              else
                 write (69, '(100f)', advance='no') save_step_pr(mm, m)
              end if
           end do
        end do
        close (69)
        open (unit=69, file = trim(fstr)//'/'//trim(chareca)//'eps.txt', action='write')
        do mm = 1, npart*nintmh
           do m = 1, npara
              if (m == npara) then
                 write (69, '(100f)') eps(mm, m)
              else
                 write (69, '(100f)', advance='no') eps(mm, m)
              end if
           end do
        end do
        close (69)

        open (unit=69, file = trim(fstr)//'/'//trim(chareca)//'randomblocks.txt', action='write')
        do mm = 1, npara
           write (69, '(i3)') ind(mm)
        end do
        close (69)

        open (unit=69, file = trim(fstr)//'/'//trim(chareca)//'break_points.txt', action='write')
        do mm = 1, nblocks+1
           write (69, '(i3)') break_points(mm)
        end do
        close (69)

     end if
     !---------------------------------------------------------------------------------
     ! End Mutation step.
     !---------------------------------------------------------------------------------

     if (rank == 0) then

        arate(i,:) = sum(acptsim,dim=2)/(1.0_wp*npart*nintmh)
        print*,''
        write(*,'(a,f5.2,f5.2,f5.2)') ' acceptance rate [min,mean,max]: ', &
             minval(arate(i,pmaskinv)), sum(arate(i,pmaskinv))/neffpara, maxval(arate(i,pmaskinv))
        write(*,'(a,f,f)'), '(unweighted) average likelihood + posterior: ', &
             sum(loglh)/(npart*1.0_wp), sum(loglh+prio)/(npart*1.0_wp)
        print*, 'mdd estimate', sum(log(zi(1:i)))

        ahat = sum(arate(i,pmaskinv))/neffpara

        ! See equation in Algorithm 3.
        scale = scale * (0.95_wp + 0.10*exp(16.0_wp*(ahat - 0.25_wp)) / (1.0_wp + &
             exp(16.0_wp*(ahat - 0.25_wp))))
        print*, 'adjusting chat', scale

        ! RECA: print every period
        !if (mod(i,10) == 0) then
        call write_files(i)
        !end if
     endif

     call mpi_bcast(scale, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
  end do

  if (rank == 0) then
     close(3)
     rfile = trim(fstr)//'/'//'resamp.txt'
     essfile = trim(fstr)//'/'//'ESS.txt'

     !------------------------------------------------------------------------------------
     ! write final output
     !------------------------------------------------------------------------------------
     open(1,file=rfile,action='write')
     open(2,file=essfile,action='write')

     do i = 1,nphi
        write(1,'(i)') resamp(i)
        write(2,'(f)') ESS(i)
     end do

     close(1)
     close(2)
     deallocate(parasim, wtsim, u, uu, arate, loglh, prio, incwt, wtsq, paraind, &
          acptsim,resamp, save_mixr, save_step_lik, save_step_pr) ! RECA
  end if

  deallocate(nodepara,nodewt,nodeloglh,nodeprio,nodeacpt, nodeeps, nodevar, &
       nodemixr, nodesteplik, nodesteppr) ! RECA
  deallocate(pmaskinv, eps, ind, break_points)
  deallocate(ESS, zi)
  call mpi_barrier(MPI_COMM_WORLD, mpierror)
  call mpi_finalize(mpierror)

contains



































  !-----------------------------------------------------------------------------
  !
  ! HELPER PROCEDURES
  !
  !-----------------------------------------------------------------------------
  subroutine transition_kernel(lk, pri, p, a, e, neval, vvar, mix, step_p0, step_lik)

    ! Does MH Transition Kernel
    integer, intent(in) :: neval
    integer, intent(inout) :: a(npara, neval)
    real(wp), intent(inout) :: lk(neval), p(npara, neval), pri(neval), e(neval*nintmh, neffpara)

    ! RECA
    real(wp), intent(inout) :: mix(nblocks*nintmh, neval)
    real(wp), intent(inout) :: step_p0(nblocks*nintmh*npara, neval)
    real(wp), intent(inout) :: step_lik(nblocks*nintmh, neval)
    ! RECA

    real(wp), intent(in):: vvar(npara, npara)
    real(wp) :: p0(npara), lik0, pr0, p1(npara), lik1, pr1
    real(wp) :: q0, q1
    real(wp) :: mixr, alp

    ! for hack
    real(wp) :: vvar2(npara,npara), sigi, zstat, ind_pdf

    real(wp), allocatable :: bvar(:,:)
    integer :: bsize
    integer, allocatable :: bi(:)
    integer :: j, k, bj, info

    a = 0
    do j = 1, neval
       p0 = p(:,j)
       lik0 = lk(j)
       pr0 = pri(j)

       vvar2 = vvar

       do bj = 1, nblocks
          bsize = break_points(bj+1) - break_points(bj)

          allocate(bi(bsize), bvar(bsize, bsize))

          bi = ind(break_points(bj)+1:break_points(bj+1))

          bvar = vvar2(bi, bi)

          call dpotrf('l', bsize, bvar, bsize, info)

          if (info .ne. 0) then
             print*,'Cholesky failed.', info
          end if


          do k = 1, nintmh
             p1 = p0

             ! RECA: SHUT THIS DOWN
             !call random_number(mixr)
             mixr = 1.0 ! 666

             ! RECA storing random number
             mix((bj-1)*(nintmh) + k, j) = mixr

             q0 = 0.0_wp
             q1 = 0.0_wp

             if (mixr < kapp) then
                ! standard RWMH
                ! recall that 'scale' has been baked into the RNG.
                p1(bi) = p1(bi) + matmul(bvar, e((j-1)*nintmh + k, bi))

             elseif (mixr < (kapp + (1.0_wp - kapp)/2.0_wp)) then
                ! diffuse RWMH with no off diagonals
                ! p1(bi) = p1(bi) + 3.0_wp*matmul(bvar, e((j-1)*nintmh + k,bi))
                ! recall that 'scale' has been baked into the RNG.
                do ii = 1,bsize
                   p1(bi(ii)) = p1(bi(ii)) + 3.0_wp*sqrt(vvar2(bi(ii),bi(ii)))*e((j-1)*nintmh+k,bi(ii))
                end do

             else
                ! independence sampler
                p1(bi) = mu(bi) + matmul(bvar, e((j-1)*nintmh + k,bi))
             endif

             ! mixture proposal densities
             q0 = kapp*exp(mvnormal_pdf(p0(bi),p1(bi),scale*bvar))
             q1 = kapp*exp(mvnormal_pdf(p1(bi),p0(bi),scale*bvar))
             ind_pdf = 1.0_wp
             do ii = 1,bsize
                sigi = sqrt(vvar2(bi(ii),bi(ii)))
                zstat = (p0(bi(ii)) - p1(bi(ii))) / sigi
                ind_pdf = ind_pdf/(sigi*sqrt(2.0_wp*3.1415))*exp(-0.5_wp*zstat**2)
                !q1 = q1 + (1.0_wp - kapp)/2.0_wp*1.0_wp/(sigi*sqrt(2.0_wp*M_PI))*exp(-0.5_wp*zstat**2)
             end do

             q0 = q0 + (1.0_wp - kapp)/2.0_wp*ind_pdf
             q1 = q1 + (1.0_wp - kapp)/2.0_wp*ind_pdf

             q0 = q0 + (1.0_wp - kapp)/2.0_wp*exp(mvnormal_pdf(p0(bi), mu(bi), scale*bvar))
             q1 = q1 + (1.0_wp - kapp)/2.0_wp*exp(mvnormal_pdf(p1(bi), mu(bi), scale*bvar))

             q0 = log(q0)
             q1 = log(q1)

             pr1 = pr(p1)
             lik1 = lik(p1)

             alp = exp(phi(i)*(lik1 - lik0) + pr1 - pr0 + q0 - q1)

!             print*, nodeu((i-2)*neval*nintmh*nblocks + (j-1)*nintmh*nblocks + (bj-1)*nintmh +k), alp, (i-2)*neval*nintmh*nblocks + (j-1)*nintmh*nblocks + (bj-1)*nintmh +k, lik1,lik0,pr1,pr0
!             print*,e((j-1)*nintmh + k,bi)
             if (nodeu((i-2)*neval*nintmh*nblocks + (j-1)*nintmh*nblocks + (bj-1)*nintmh +k) < alp) then
                p0 = p1
                lik0 = lik1
                pr0 = pr1
                a(bi,j) = a(bi,j) + 1
             endif
             ! RECA
             step_p0((bj-1)*(nintmh)*npara + (k-1)*npara + 1:(bj-1)*nintmh*npara + k*npara, j) = p0
             step_lik((bj-1)*(nintmh)*npara + k, j) = lik0

          end do

          p(:,j) = p0
          lk(j) = lik0
          pri(j) = pr0
          deallocate(bi, bvar)
       end do
    end do
  end subroutine transition_kernel


  subroutine initial_draw_from_prior(lk, pri, p, wt, neval, npara)
    ! Computes the initial draw from the prior.

    integer, intent(in) :: neval, npara
    real(wp), intent(inout) :: lk(neval), p(npara, neval), wt(neval), pri(neval)

    integer :: i
    real(wp) :: p2(2,npara)

    real(wp) :: lik0

    do i = 1, neval
       lik0 = lik(p(:,i))
       !print*,lik0
       if (isnan(lik0)) then
          lik0 = -10.0_wp**11
       endif

       do while (lik0 < -10.0_wp**9)

          p2 = priorrand(2, pshape, pmean, pstdd, pmask, pfix)

          lik0 = lik(p2(1,:))

          if (isnan(lik0)) then
             lik0 = -10000000000000.0_wp
          endif

          p(:,i)= p2(1,:)
          wt(i) = pr(p(:,i))


       end do

       lk(i) = lik0
       pri(i) = pr(p(:,i))

       if (mod(i, 100) == 0) then
          write(*,'(a$)') '.'
       endif
    end do

  end subroutine initial_draw_from_prior




  subroutine scatter(everything)
    ! spread data to nodes

    logical, intent(in), optional :: everything
    logical :: everything_0

    if (not(present(everything))) then
       everything_0 = .false.
    else
       everything_0 = everything
    endif

    ! now we are ready to start the SMC algorithm
    call mpi_scatter(parasim, ngap*npara, MPI_DOUBLE_PRECISION, nodepara, ngap*npara, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
    call mpi_scatter(wtsim, ngap, MPI_DOUBLE_PRECISION, nodewt, ngap, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
    call mpi_scatter(loglh, ngap, MPI_DOUBLE_PRECISION, nodeloglh, ngap, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
    call mpi_scatter(prio, ngap, MPI_DOUBLE_PRECISION, nodeprio, ngap, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)

    if (everything_0 == .true.) then

       ! RECA
       call mpi_scatter(save_mixr, ngap*nintmh*nblocks, MPI_DOUBLE_PRECISION, nodemixr, ngap*nintmh*nblocks, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
       call mpi_scatter(save_step_lik, ngap*nintmh*nblocks, MPI_DOUBLE_PRECISION, nodesteplik, ngap*nintmh*nblocks, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
       call mpi_scatter(save_step_pr, ngap*nintmh*nblocks*npara, MPI_DOUBLE_PRECISION, nodesteppr, ngap*nintmh*nblocks*npara, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
       !!!!

       call mpi_scatter(acptsim, ngap*npara, MPI_INTEGER, nodeacpt, ngap*npara, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)
       call mpi_scatter(eps, ngap*nintmh*npara, MPI_DOUBLE_PRECISION, nodeeps, ngap*nintmh*npara, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
       call mpi_bcast(var, npara*npara, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
       call mpi_bcast(mu, npara, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
       call mpi_bcast(ind, neffpara, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)
       call mpi_bcast(break_points, nblocks+1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)

    end if
  end subroutine scatter




  subroutine gather(everything)
    ! send data to master
    logical, intent(in), optional :: everything
    logical :: everything_0
    if (not(present(everything))) then
       everything_0 = .false.
    else
       everything_0 = everything
    endif

    call mpi_gather(nodewt, ngap, MPI_DOUBLE_PRECISION, wtsim, ngap, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )
    call mpi_gather(nodeloglh, ngap, MPI_DOUBLE_PRECISION, loglh, ngap, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )
    call mpi_gather(nodeprio, ngap, MPI_DOUBLE_PRECISION, prio, ngap, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )
    call mpi_gather(nodepara, ngap*npara, MPI_DOUBLE_PRECISION, parasim, ngap*npara, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )

    ! RECA
    call mpi_gather(nodemixr, ngap*nintmh*nblocks, MPI_DOUBLE_PRECISION, save_mixr, ngap*nintmh*nblocks, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )
    call mpi_gather(nodesteplik, ngap*nintmh*nblocks, MPI_DOUBLE_PRECISION, save_step_lik, ngap*nintmh*nblocks, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )
    call mpi_gather(nodesteppr, ngap*nintmh*nblocks*npara, MPI_DOUBLE_PRECISION, save_step_pr, ngap*nintmh*nblocks*npara, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )


    if (everything_0 == .true.) then
       call mpi_gather(nodeacpt, ngap*npara, MPI_INTEGER, acptsim, ngap*npara, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror )
    endif

  end subroutine gather



  subroutine write_hyperparameters(ii)
    !! save hyperparameters mu, var, c

    integer, intent(in) :: ii

    character(len=244) :: chari, blockfile,varfile,meanfile,scalefile

    integer :: ij


    write(chari,'(i3.3)') ii

    blockfile = trim(fstr)//'/'//trim(chari)//'blocksim.txt'
    varfile = trim(fstr)//'/'//trim(chari)//'var.txt'
    meanfile = trim(fstr)//'/'//trim(chari)//'mean.txt'
    scalefile = trim(fstr)//'/'//trim(chari)//'scale.txt'


    open(1, file=blockfile, action='write')
    open(2, file=varfile, action='write')
    open(4, file=meanfile, action='write')
    open(5, file=scalefile, action='write')


    write(5, '(1f)') scale

    do ij = 1,npara
       write(4,'(100f)') mu(ij)
       write(2,'(100f)') var(ij,:)
    end do

    write(1,'(100i)') ind
    write(1,'(100i)') break_points

    close(1)
    close(2)
    close(4)
    close(5)

  end subroutine write_hyperparameters

  subroutine load_hyperparameters(ii)

    integer, intent(in) :: ii

    integer :: ij

    character(len=244) :: chari, blockfile,varfile,meanfile,scalefile

    write(chari,'(i3.3)') ii

    blockfile = trim(hyp_fstr)//'/'//trim(chari)//'blocksim.txt'
    varfile = trim(hyp_fstr)//'/'//trim(chari)//'var.txt'
    meanfile = trim(hyp_fstr)//'/'//trim(chari)//'mean.txt'
    scalefile = trim(hyp_fstr)//'/'//trim(chari)//'scale.txt'


    open(1, file=blockfile, action='read')
    open(2, file=varfile, action='read')
    open(4, file=meanfile, action='read')
    open(5, file=scalefile, action='read')

    read(5,*) scale

    do ij=1,npara
       read(4,*) mu(ij)
       read(2,*) var(ij,:)
    end do

    read(1,*) ind
    read(1,*) break_points

    close(1)
    close(2)
    close(4)
    close(5)

  end subroutine load_hyperparameters

  subroutine write_files(ii)
    ! write files
    integer, intent(in) :: ii

    character(len=144) :: chari

    write(chari,'(i3.3)') ii

    parafile = trim(fstr)//'/'//trim(chari)//'parasim.txt'
    wghtfile = trim(fstr)//'/'//trim(chari)//'wtsim.txt'
    llikfile = trim(fstr)//'/'//trim(chari)//'liksim.txt'
    postfile = trim(fstr)//'/'//trim(chari)//'postsim.txt'

    open(1, file=parafile, action='write')
    open(2, file=wghtfile, action='write')
    open(4, file=llikfile, action='write')
    open(5, file=postfile, action='write')
    do j = 1, npart

       write(1,'(100f)') parasim(:,j)
       write(2,'(1f20.7)') wtsim(j)
       write(4,'(1f20.7)') loglh(j)
       write(5,'(1f20.7)') loglh(j)+prio(j)

    end do

    close(1)
    close(2)
    close(4)
    close(5)

    write(*,*) 'Wrote files '//trim(chari)//'{parasim,wtsim,liksim,postsim}.txt'
    print*,''
  end subroutine write_files

  subroutine read_in_from_files()
    ! read in prior, data, translation from textfiles
    open(1, file=priorfile, status='old', action='read')
    do i = 1, npara
       read(1, *) pshape(i), pmean(i), pstdd(i), pmask(i), pfix(i)
    end do
    close(1)

    open(1, file=datafile, status='old', action='read')
    do i = 1, nobs
       read(1, *) YY(:,i)
    end do
    close(1)

    open(1, file=transfile, status='old', action='read')
    do i = 1, npara
       read(1, *) trspec(:,i)
    end do
    close(1)
  end subroutine read_in_from_files

  subroutine print_help()
    ! help file
    print '(a)', 'smc_driver_mpi -- sequential monte carlo estimation for DSGE models'
    print '(a)', '      by Ed Herbst [edward.p.herbst@frb.gov]'
    print '(a)', ''
    print '(16a,a)', 'current model  : ', mname
    print '(18a,a)', 'current basedir: ', basedir
    print '(a)', ''
    print '(a)', 'usage: mpirun -n NPROC ./smc_driver_mpi [OPTIONS]'
    print '(a)', ''
    print '(a)', 'options:'
    print '(a)', ''
    print '(a)', '-p, --nphi [N]        sets nphi = N                       DEFAULT = 500'
    print '(a)', '-n, --npart [N]       sets nparticle = N                  DEFAULT = 4096'
    print '(a)', '-m, --nintmh [N]      sets the number of mh steps = N     DEFAULT = 1'
    print '(a)', '-b, --bend [LAM]      bends phi^LAM                       DEFAULT = 1'
    print '(a)', '-u,                   does not use a mixture proposal     DEFAULT not selected'
    print '(a)', '-g, --geweke          does a partial posterior SMC        DEFAULT not selected'
    print '(a)', '-o, --nblocks [N]     sets number of blocks = N           DEFAULT = 1'
    print '(a)', '-i, --initrep [N]     sets the trial number = N           DEFAULT = 1'
    print '(a)', ''
  end subroutine print_help
end program smc_driver
