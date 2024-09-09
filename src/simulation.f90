!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_uns_class,   only: hypre_uns
   use hypre_str_class,   only: hypre_str
   use incomp_class,      only: incomp
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use rotorDisk_class,      only: rotorDisk
   use blade_class,          only: blade
   implicit none
   private

   !> Single-phase incompressible flow solver and corresponding time tracker
   type(incomp),      public :: fs
   type(timetracker), public :: time
   type(hypre_uns),     public :: ps
   type(hypre_str),     public :: vs
   type(rotorDisk),      public :: rd
   type(blade),          public :: bl

   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,rdfile

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: rho

   !> Fluid viscosity
   real(WP) :: visc

   real(WP) :: U1, U2, U3, U4
   real(WP) :: P1, P2, P3, P4


contains

   !> Function that localizes the left (x-) of the domain
   function left_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function left_of_domain

   !> Function that localizes the right (x+) of the domain
   function right_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_of_domain

   function bottom_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function bottom_of_domain

   function top_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function top_of_domain

   function front_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function front_of_domain

   function back_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function back_of_domain

   subroutine getSpecificVelocity()
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM, MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i,j,k
      real(WP) :: U1_loc, U2_loc, U3_loc, U4_loc
      real(WP) :: P1_loc, P2_loc, P3_loc, P4_loc
      integer :: ierr
      integer :: mycount1, mycount2, mycount3, mycount4
      integer :: count1, count2, count3, count4

      integer :: rotordiskIndex
      integer,  dimension(3) :: ind

      U1_loc = 0.0_WP ; U2_loc = 0.0_WP ; U3_loc = 0.0_WP ; U4_loc = 0.0_WP
      P1_loc = 0.0_WP ; P2_loc = 0.0_WP ; P3_loc = 0.0_WP ; P4_loc = 0.0_WP

      mycount1 = 0 ; mycount2 = 0 ; mycount3 = 0 ; mycount4 = 0

      ind = cfg%get_ijk_global(rd%center,[0,0,0])
      rotordiskIndex = ind(1)

      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_

               ! U1, ahead of the rotor disk
               if ((cfg%xm(i).lt.(rd%center(1)-1.5_WP+cfg%dx(i))) .and. (cfg%xm(i).gt.(rd%center(1)-1.5_WP-cfg%dx(i)))) then
                  if ((cfg%ym(j).lt.(rd%center(2)+2.5_WP*rd%maxR)) .and. (cfg%ym(j).gt.(rd%center(2)-2.5_WP*rd%maxR))) then
                     U1_loc = U1_loc + Ui(i,j,k)
                     P1_loc = P1_loc + fs%P(i,j,k)
                     mycount1 = mycount1 + 1
                  end if
               end if

               ! U4, behind the rotor disk
               if ((cfg%xm(i).lt.(rd%center(1)+2.0_WP+cfg%dx(i))) .and. (cfg%xm(i).gt.(rd%center(1)+2.0_WP-cfg%dx(i)))) then
                  if ((cfg%ym(j).lt.(rd%center(2)+4.0_WP*rd%maxR)) .and. (cfg%ym(j).gt.(rd%center(2)-4.0_WP*rd%maxR))) then
                     U4_loc = U4_loc + Ui(i,j,k)
                     P4_loc = P4_loc + fs%P(i,j,k)
                     mycount4 = mycount4 + 1
                  end if
               end if   

            end do
         end do
      end do

      call MPI_ALLREDUCE(U1_loc, U1, 1, MPI_REAL_WP, MPI_SUM, cfg%comm, ierr)
      call MPI_ALLREDUCE(U4_loc, U4, 1, MPI_REAL_WP, MPI_SUM, cfg%comm, ierr)

      call MPI_ALLREDUCE(P1_loc, P1, 1, MPI_REAL_WP, MPI_SUM, cfg%comm, ierr)
      call MPI_ALLREDUCE(P4_loc, P4, 1, MPI_REAL_WP, MPI_SUM, cfg%comm, ierr)

      call MPI_ALLREDUCE(mycount1, count1, 1, MPI_INTEGER, MPI_SUM, cfg%comm, ierr)
      call MPI_ALLREDUCE(mycount4, count4, 1, MPI_INTEGER, MPI_SUM, cfg%comm, ierr)

      U1 = U1 / real(count1, WP)
      U4 = U4 / real(count4, WP)

      P1 = P1 / real(count1, WP)
      P4 = P4 / real(count4, WP)





   end subroutine getSpecificVelocity


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         allocate(rho(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))


      end block allocate_work_arrays


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=1
      end block initialize_timetracker


      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_uns_class, only: gmres_amg
         use hypre_str_class, only: pcg_pfmg
         use incomp_class,      only: clipped_neumann, dirichlet
         use incomp_class, only: bcond
         integer :: i,j,k,n
         type(bcond), pointer :: mybc

         ! Create flow solver
         fs=incomp(cfg=cfg,name='NS solver')

         call fs%add_bcond(name='inflow',type=dirichlet,locator=left_of_domain,face='x',dir=-1,canCorrect=.false. )
         call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true. )
         call fs%add_bcond(name='bottom',type=clipped_neumann,locator=bottom_of_domain,face='y',dir=-1,canCorrect=.true. )
         call fs%add_bcond(name='top',type=clipped_neumann,locator=top_of_domain,face='y',dir=+1,canCorrect=.true. )
         ! call fs%add_bcond(name='front',type=clipped_neumann,locator=front_of_domain,face='z',dir=-1,canCorrect=.true. )
         ! call fs%add_bcond(name='back',type=clipped_neumann,locator=back_of_domain,face='z',dir=+1,canCorrect=.true. )

         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign constant density
         call param_read('Density',fs%rho)
         ! Configure pressure solver
         ps=hypre_uns(cfg=cfg,name='Pressure',method=gmres_amg,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)



         fs%U=1.0_WP ; fs%V=0.0_WP ; fs%W=0.0_WP
         ! Apply all other boundary conditions 
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)   = 1.0_WP
         end do
         call fs%apply_bcond(time%t,time%dt)

         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
      end block create_and_initialize_flow_solver

      ! Create a rotor disk
      create_and_initialize_rotorDisk: block
         use param, only: param_read

         integer :: Nr, Na          ! Number of radial and azimuthal points

         call param_read('Number of radical points',Nr)
         call param_read('Number of azimuthal points',Na)

         ! create blade object
         bl=blade(Nr=Nr,Na=Na)

         ! read blade parameters
         call param_read('Blade radius',bl%radius)
         call param_read('Blade twists',bl%twist)
         call param_read('Blade Chords',bl%chord)
         call param_read('Blade AoA',bl%aoa)
         call param_read('Blade Cd',bl%cd)
         call param_read('Blade Cl',bl%cl)

         rd=rotorDisk(cfg=cfg,bl=bl)

         call param_read('Rotor Disk number of blades',rd%nblades)
         call param_read('Rotor Disk min radius',rd%minR)
         call param_read('Rotor Disk max radius',rd%maxR)
         call param_read('Rotor Disk center',rd%center)
         call param_read('Rotor Disk axis',rd%axis)
         call param_read('Rotor Disk reference direction',rd%ref_dir)
         call param_read('Rotor Disk angular velocity',rd%omega)

         ! prepare the rotor disk
         call rd%prepareRotorDisk()

         ! call rd%setTipEffect(0.86_WP)

         rho=fs%rho
         call cfg%sync(rho)
         call rd%calculateForce(rho,Ui,Vi,Wi)    ! Get volumetric force


      end block create_and_initialize_rotorDisk


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='channel')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('viscosity',fs%visc)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('div',fs%div)
         call ens_out%add_scalar('diskArea',rd%area)
         call ens_out%add_vector('diskForce',rd%forceX,rd%forceY,rd%forceZ)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create forcing monitor
         rdfile=monitor(fs%cfg%amRoot,'rotorDisk')
         call rdfile%add_column(time%n,'Timestep number')
         call rdfile%add_column(time%t,'Time')
         call rdfile%add_column(rd%Thrust,'Thrust')
         call rdfile%add_column(rd%Torque,'Torque')
         call rdfile%add_column(rd%Power,'Power')
         call rdfile%add_column(U1,'U1')
         call rdfile%add_column(U2,'U2')
         call rdfile%add_column(U3,'U3')
         call rdfile%add_column(U4,'U4')
         call rdfile%add_column(P1,'P1')
         call rdfile%add_column(P2,'P2')
         call rdfile%add_column(P3,'P3')
         call rdfile%add_column(P4,'P4')
         call rdfile%write()
         

      end block create_monitor



   end subroutine simulation_init


   !> Time integrate our problem
   subroutine simulation_run
      use incomp_class, only: bcond
      implicit none
      integer :: i,j,k,n
      type(bcond), pointer :: mybc
      integer, dimension(3) :: index

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Add rotor disk momentum source terms
            ! call rd%calculateForce(rho,Ui,Vi,Wi)    ! Get volumetric force
            ! resU=resU+rd%forceX
            ! resV=resV+rd%forceY
            ! resW=resW+rd%forceZ
            
            ! Explicit Momentum source term
            index = cfg%get_ijk_global(rd%center,[0,0,0])
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_

                     if (i == index(1)) then
                        if ((cfg%ym(j).le.(rd%center(2)+rd%maxR)) .and. (cfg%ym(j).ge.(rd%center(2)-rd%maxR))) then
                              resU(i,j,k) = resU(i,j,k) + 100.0_WP
                        end if
                     end if
                     
                  end do
               end do
            end do


            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW


            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply other boundary conditions on the resulting fields
            call fs%get_bcond('inflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               fs%U(i,j,k)   = 1.0_WP
            end do
            call fs%apply_bcond(time%t,time%dt)

            ! Solve Poisson equation
            call fs%correct_mfr()
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho
            fs%V=fs%V-time%dt*resV/fs%rho
            fs%W=fs%W-time%dt*resW/fs%rho

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         call getSpecificVelocity()

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()

         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         call rdfile%write()  

      end do

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker

      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)

   end subroutine simulation_final





end module simulation
