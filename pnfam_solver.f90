!------------------------------------------------------------------------------
!> This module contains the heart of the pnfam solver. It calls the setup routines
!> runs the iterative solver, and writes the output to file.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module pnfam_solver
   use pnfam_setup
   use type_extfield
   use type_bigblockmatrix

   implicit none

   !----------------------------------QQst      modify18
   character(100), public,  parameter :: qq_io_path = "/21dayscratch/scr/q/u/qunqun/pn_io_1/"
   !----------------------------------QQend     modify18



   integer, parameter, private :: dp = kind(1d0)

   logical, private :: use_diagonal_blocks, bminus

   ! External field
   type(bigblockmatrix), private :: Fsp, Fqp
   type(bigblockmatrix), dimension(:), allocatable, private :: Gsp, Gqp
   ! HFB solution
   type(bigblockmatrix), private :: Wn, Wp
   ! FAM matrices
   type(bigblockmatrix), private :: dHqp_re, dHqp_im, dRqp_re, dRqp_im
   type(bigblockmatrix), private :: dHsp_re, dHsp_im, dRsp_re, dRsp_im
   type(bigblockmatrix), private :: Greens_re, Greens_im, T

contains
   !---------------------------------------------------------------------------
   ! Subroutine for solving the FAM equations for one value of energy
   !---------------------------------------------------------------------------
   subroutine pnfam_solve
      use pnfam_broyden
      use pnfam_txtoutput
      use pnfam_storage, only : fam_io


      implicit none
      integer :: iexit=1, ierr=0

      ! Initialize FAM inputs (as blockmatrix types, where applicable)
      ! - Logging procedures
      ! - HFB solution (E,U,V)
      ! - External field (F in SP basis)
      ! - Routine for calculating dH
      !---------------------------------------------------------------
      call setup_pnfam
      call openlog(output_txtfile, print_stdout)
      call print_main_header
      call setup_extfield

      ! Prepare ifam
      !---------------------------------------------------------------
      ! 1. Allocate Block Matices:
      !    - F_QP, G_QP, dH_QP(H20,H02,H11,H11t), dH_SP(h,d),
      !      dR_QP(X,Y,P,Q), dR_SP(rho,kap), Greens(X,Y,P,Q), T(1ff,ff)
      ! 2. Transform F to QP basis
      ! 3. Define block structures (mostly based on F)
      ! 4. Populate statmatrix T and Greens function
      !---------------------------------------------------------------
      call init_pnfam_solver

      ! Initialize fam amplitudes
      !---------------------------------------------------------------
      call fam_io(1,iexit)
      if (iexit > 0) then
         si = 1
         call set_val_bbm(dRqp_re,0d0)
         call set_val_bbm(dRqp_im,0d0)
         call set_broinout_from_dRqp(dRqp_re, dRqp_im,'i')
      else
         call set_dRqp_from_broin(dRqp_re, dRqp_im)
      end if

      ! Run the iterative solver
      !---------------------------------------------------------------
      call log_header
      call ifam

      ! Write the result to footer and store result
      !---------------------------------------------------------------
      call log_result
      call fam_io(2,iexit)
      call closelog

      !call transition_density

   end subroutine pnfam_solve

   !---------------------------------------------------------------------------
   ! Subroutine for solving the FAM equations for one value of energy. Before
   ! calling this the first time, the module data structures must be initialized
   ! with init_pnfam_solver.
   !---------------------------------------------------------------------------
   subroutine ifam
      use pnfam_constants, only : iu, pi, get_timer
      use pnfam_broyden

     !---------------------------QQst
     ! qqphonon
     !-----------------
      use pnfam_WX
     !---------------------------Qend



      implicit none

      ! Finite temperature
      real(dp)    :: ft_real_factor
      complex(dp) :: ft_cmplx_factor

      integer :: iter, ixterm


      real(dp) :: re_s, im_s
      real(dp) :: time1, time2, timei, timei0


      !--------------------------------------------QQst
      ! variable for assigning X,Y to matX and matY
      !---------------------------------------
     !integer, parameter :: dp = kind(1.0d0)
      integer, parameter :: dc = kind(cmplx(1.0,1.0,kind=dp))
      integer, allocatable :: Ph_readin_select(:,:)
      integer, allocatable :: U_select(:,:),V_select(:,:)
      integer, allocatable :: X_select(:,:),Y_select(:,:)
      integer, allocatable :: Ph_db(:)
      integer :: Ph_nb, Ph_NumberofX
      integer :: Ph_sum_db
      integer, allocatable :: Ph_sum_db_colm1(:)          !!!Q: dimension is iHFB_NB, record sum of nd of first column minus 1 blocks.
      Complex(dc), allocatable  :: Ph_matX(:,:),Ph_matY(:,:)
      Complex(dc), allocatable  :: Ph_X(:), Ph_Y(:)       !!!Q: maybe not necessary.
      !---------------------------------------
      ! variable for assigning matN11 and matP11
      !---------------------------------------
      !
      integer, allocatable :: WX_readin_select(:,:)
      integer, allocatable :: WX_id(:)
      integer :: WX_iHFB_NB, WX_NumberofPH11,WX_NumberofPH11_maxval
      integer :: WX_sum_nd
      integer, allocatable :: WX_sum_nd_colm1(:)          !!!Q: dimension is iHFB_NB, record sum of nd of first column minus 1 blocks.
      real(dp) :: WX_realomega 
      real(dp) :: W_omega 
      !---------------------------------------
      Complex(Kind(1.000D0)), allocatable  :: WX_NH11_E101(:),WX_PH11_E101(:)
      Complex(Kind(1.000D0)), allocatable  :: WX_BNH11_E101(:),WX_BPH11_E101(:)  !!!Q: bmod01
      !
      Complex(Kind(1.000D0)), allocatable  :: WX_matN11_E101(:,:),WX_matP11_E101(:,:)
      Complex(Kind(1.000D0)), allocatable  :: WX_matBN11_E101(:,:),WX_matBP11_E101(:,:)  !!!Q: bmod01
      !
      Complex, allocatable  :: WX_matpnbasisN11_E101(:,:),WX_matpnbasisP11_E101(:,:) !!!Q:bmod01    
     !real(dp), allocatable  :: WX_matpnbasisN11_E101(:,:),WX_matpnbasisP11_E101(:,:)
      integer, allocatable :: WX_pnbasis_select(:,:)   !modifyX01
      !---------------------------------------
 
      ! variable used by both matX, matY and matN11, matP11
      !---------------------------------------
      integer :: sum_block
      integer :: qqj, qqi,qqk, ct1, ct2, ct3      !!!Q: ct means count, for any
      integer :: inbcc                        !!!Q: inside block column count, used to decide which column to be assigned.
      integer :: inbrc                        !!!Q: inside block row count, used to decide whcih row element to be assigned.
      integer :: inbrr                        !!!Q: used to record how many elements already assigned to a row inside a block.
      logical :: Ph_test=.True.               !!!Q: logical to open or close phonon calculation.
      complex,parameter :: cunit=cmplx(0.0_dp,1.0_dp,kind=dp)
      complex,parameter :: runit=cmplx(1.0_dp,0.0_dp,kind=dp)
      integer :: nph, nph_rec   !!! number of phonons 

     !Complex(Kind(1.000D0)), allocatable  :: WX_NH11_store(:,:),WX_PH11_store(:,:)
     !Complex(Kind(1.000D0)), allocatable  :: WX_BNH11_store(:,:),WX_BPH11_store(:,:)  !!!Q: bmod01
      real(dp), allocatable :: WX_realomega_store(:) 
      integer, allocatable :: WX_iHFB_NB_store(:), WX_NumberofPH11_store(:)
      integer, allocatable :: WX_id_store(:,:)
      integer, allocatable :: WX_readin_select_store(:,:,:)
      integer, allocatable :: iso_store(:), L_store(:), k_store(:)



      Complex(Kind(1.000D0)), allocatable  :: WX_NH11_store1D(:)
      Complex(Kind(1.000D0)), allocatable  :: WX_BNH11_store1D(:)
      Complex(Kind(1.000D0)), allocatable  :: WX_PH11_store1D(:)
      Complex(Kind(1.000D0)), allocatable  :: WX_BPH11_store1D(:)


      !---------------------------------------------
      ! matrix PX, PPX,  
      !---------------------------------------------
      !interm variable matrix, used for both WX and WY
   !  Complex(Kind(1.000D0)), allocatable  :: Ph_matPdF(:,:),Ph_matPdFg(:,:), Ph_matFNc(:,:),Ph_matFNcg(:,:)
      ! matrix distinguished.
      !---------------------------------------------
   !  Complex(Kind(1.000D0)), allocatable  :: Ph_matPPdXg(:,:), Ph_matXNcgN(:,:), Ph_matPXNcg(:,:), Ph_matPdXgN(:,:)
   !  Complex(Kind(1.000D0)), allocatable  :: Ph_matPPdYg(:,:), Ph_matYNcgN(:,:), Ph_matPYNcg(:,:), Ph_matPdYgN(:,:)
   !  Complex(Kind(1.000D0)), allocatable  :: Ph_matWX(:,:), Ph_matWY(:,:)             !modify32
      !---------------------------------------------
      real, allocatable :: Ph_ReWX(:),Ph_ImWX(:)  
      real, allocatable :: Ph_ReWY(:),Ph_ImWY(:)  
      !---------------------------------------------
      Complex, allocatable  :: WgreenX(:), WgreenY(:)                                             
      !---------------------------------------------
 
      real :: Ph_mult=1.0
      real :: Delta = +0.1
      !-------------------------QQst    modify30
      real :: si_mini=100.0
      integer :: iter_start, iter_end
      integer, parameter :: iter_endmax = 100
      !-------------------------QQend   modify30
      complex, dimension(1) :: cstr_mini
      !--------------------------
      integer :: which_iso
      integer :: which_L
      integer :: which_k
      integer :: which_ph
      character(len=10)  :: iso_id
      character(len=10)  :: L_id
      character(len=10)  :: k_id
      character(len=10)  :: ph_id
      logical :: file_exists
     
      open(48,file='Ph_converge',access='append')    ! modify18
     !open(49,file='Ph_from_simplex_to_pnfam_dH11')
      open(50,file='50_checke')    ! modify18
     !
      !--------------------------------------------QQend

      !---------------------------------------------------------------QQst  modify20
      ! allocate for Phonon calculation.
      !---------------------
      if(Ph_test) then			!!!Q: prepare before iter loop
          !--------------------------------
	  !--------------------------------
	  ! readin binfo from exterfield, assign array X and Y to matX and matY
	  !--------------------------------
	  open(300,file='Ph_binary_output_extfield_GT0',status='unknown',form='unformatted')
	  Read(300) Ph_nb
	  Read(300) Ph_NumberofX                    !!!Q: read
	  !
	  if(allocated(Ph_db)) deallocate(Ph_sum_db_colm1,Ph_db)
	  Allocate(Ph_sum_db_colm1(Ph_nb),Ph_db(Ph_nb))
	  ! 
	  Read(300) Ph_db                 !!!Q: pnfam_solver.f90 does't have db(:)
	  Read(300) Ph_sum_db
	  Read(300) Ph_sum_db_colm1
	  !
	  !
	  if(allocated(Ph_readin_select)) deallocate(Ph_readin_select)
	  Allocate(Ph_readin_select(Ph_nb,Ph_nb))
	  !
	  Read(300) Ph_readin_select                   !!!Q: read
	  close(300)

	  Allocate(U_select(Ph_nb,Ph_nb), V_select(Ph_nb,Ph_nb), &
                   X_select(Ph_nb,Ph_nb), Y_select(Ph_nb,Ph_nb))

         U_select = 0
         V_select = 0
         X_select = 0
         Y_select = 0
         do qqi = 1, Ph_nb  
             do qqj = 1, Ph_nb
                 if(qqi.eq.(qqj-Ph_nb/2)) then
                     V_select(qqi,qqj) = 1
                 end if
                 if((qqi-Ph_nb/2).eq.qqj) then
                     V_select(qqi,qqj) = 1
                 end if
             end do 
         end do
         !
         X_select= matmul(Ph_readin_select, V_select)
         Y_select= matmul(V_select, Ph_readin_select)
         open(20,file='X_select')
         do qqi = 1, Ph_nb
            do qqj = 1, Ph_nb
              write(20,*)qqi, qqj, U_select(qqi,qqj), V_select(qqi,qqj), &
                        X_select(qqi,qqj)
            end do
         end do
         close(20)


          !-----------------------
          ! Wgreens         
          !----------------------
	  if(allocated(WgreenX)) deallocate(WgreenX, WgreenY)
	  Allocate(WgreenX(Ph_sum_db**2), WgreenY(Ph_sum_db**2))
          WgreenX = 0
          WgreenY = 0
	  !
         !Call WgreenX_func(Ph_sum_db, Ep, En, omega, WX_realomega, width, WgreenX) 
         !Call WgreenY_func(Ph_sum_db, Ep, En, omega, WX_realomega, width, WgreenY) 
          !
          Allocate(Ph_ReWX(Ph_NumberofX), Ph_ImWX(Ph_NumberofX))
          Ph_ReWX=0
          Ph_ImWX=0
 
          Allocate(Ph_ReWY(Ph_NumberofX), Ph_ImWY(Ph_NumberofX))
          Ph_ReWY=0
          Ph_ImWY=0
 
	  !----------------------------
	  ! preserve below to check
	  !----------------------------
	 !sum_block=0
	 !do qqj=1,Ph_nb
	 !   do qqi=1,Ph_nb
	!if(Ph_readin_select(qqj,qqi)==1) then 
	!   sum_block=sum_block+Ph_db(qqj)*Ph_db(qqi)             !!!Q: number of elements sum block 
	!  !write(34,*) qqj, qqi, Ph_db(qqj), Ph_db(qqi), sum_block
	!end if
	!    end do
	! end do
          !
          !-----------------------------
          ! allocate Ph_matX and Ph_matY 
          !-----------------------------
          Allocate(Ph_matX(Ph_sum_db,Ph_sum_db),Ph_matY(Ph_sum_db,Ph_sum_db))
          Ph_matX=0
          Ph_matY=0



         nph = 0
         !=================================================

         !=================================================
         ! find out how many phonons, and allocate storage memory.
         !=================================================
        do which_iso=  0,1
           write(iso_id,'(i0)') which_iso
             do which_L = 0,7
               write(L_id,'(i0)') which_L
               do which_k = 0, which_L
                 do which_ph = 1,15
                   write(k_id,'(i0)') which_k
                   write(ph_id,'(i0)') which_ph
                   !if(which_k.eq.0) then
                   !   Ph_mult = 1                   !!!Q: apply for deformed nuclear.          
                   !else
                   !   Ph_mult = 2  
                   !end if
                 
                

	  inquire(file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),exist=file_exists)    !modify18
      if(file_exists) then
         nph = nph + 1
      endif
    
                 end do
               end do
             end do
        end do

        nph_rec = nph

        write(50, *) 'nph:', nph


	  if(allocated(WX_realomega_store)) deallocate(WX_realomega_store, WX_iHFB_NB_store, WX_NumberofPH11_store)
	  allocate(WX_realomega_store(nph), WX_iHFB_NB_store(nph), WX_NumberofPH11_store(nph))
      if(allocated(iso_store)) deallocate(iso_store, L_store, k_store)
      allocate(iso_store(nph), L_store(nph), k_store(nph))

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        nph = 0
        do which_iso=  0,1
        write(iso_id,'(i0)') which_iso
        do which_L = 0,7
        write(L_id,'(i0)') which_L
        do which_k = 0, which_L
           do which_ph = 1,15
              write(k_id,'(i0)') which_k
              write(ph_id,'(i0)') which_ph
                if(which_k.eq.0) then
                   Ph_mult = 1                   !!!Q: apply for deformed nuclear.          
                else
                   Ph_mult = 2  
                end if
                

	  inquire(file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),exist=file_exists)    !modify18
      if(file_exists) then
      nph = nph  + 1
	  open(110,file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),status='unknown',form='unformatted')       !modify18
          Read(110) WX_realomega
	  Read(110) WX_iHFB_NB                         !!!Q: read


	  if(allocated(WX_id)) deallocate(WX_id, WX_sum_nd_colm1)
	  allocate(WX_id(WX_iHFB_NB), WX_sum_nd_colm1(WX_iHFB_NB))

          if(allocated(WX_readin_select))  deallocate(WX_readin_select)
	  Allocate(WX_readin_select(WX_iHFB_NB,WX_iHFB_NB))


	  Read(110) WX_NumberofPH11                    !!!Q: read
	  Read(110) WX_id                              !!!Q: read
	  Read(110) WX_readin_select                   !!!Q: read
	  !
	  close(110)
 
      WX_realomega_store(nph) = WX_realomega
      WX_iHFB_NB_store(nph) = WX_iHFB_NB
      WX_NumberofPH11_store(nph) = WX_NumberofPH11
 
      iso_store(nph) =  which_iso
      L_store(nph) =  which_L 
      k_store(nph) =  which_k 

      end if    ! file_exists
         end do ! which_ph
         end do ! which_k
         end do ! which_L
         end do ! which_iso

    WX_NumberofPH11_maxval = maxval(WX_NumberofPH11_store)
    write(50,*) 'WX_realomega_store' ,WX_realomega_store
    write(50,*) 'WX_iHFB_NB_store' ,WX_iHFB_NB_store 
    write(50,*) 'WX_NumberofPH11_store',WX_NumberofPH11_store
    write(50,*) 'WX_NumberofPH11_maxval',WX_NumberofPH11_maxval


	  if(allocated(WX_id_store)) deallocate(WX_id_store, WX_readin_select_store)
	  allocate(WX_id_store(nph, WX_iHFB_NB_store(1)), WX_readin_select_store(nph, WX_iHFB_NB_store(1),WX_iHFB_NB_store(1)))


	 !if(allocated(WX_NH11_store)) deallocate(WX_NH11_store, WX_BNH11_store, WX_PH11_store, WX_BPH11_store)
	 !allocate(WX_NH11_store(nph, WX_NumberofPH11_maxval), WX_BNH11_store(nph, WX_NumberofPH11_maxval), & 
     !          WX_PH11_store(nph, WX_NumberofPH11_maxval), WX_BPH11_store(nph, WX_NumberofPH11_maxval))

      if(allocated(WX_NH11_store1D)) deallocate(WX_NH11_store1D,WX_BNH11_store1D,WX_PH11_store1D,WX_BPH11_store1D)
      allocate(WX_NH11_store1D(nph*WX_NumberofPH11_maxval))
     !allocate(WX_BNH11_store1D(nph*WX_NumberofPH11_maxval))
      allocate(WX_PH11_store1D(nph*WX_NumberofPH11_maxval))
     !allocate(WX_BPH11_store1D(nph*WX_NumberofPH11_maxval))

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! loop phonon binfo once more to populate WX_id_store,
    ! WX_readin_select_store,  WX_?H11_store
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        nph = 0
        do which_iso=  0,1
        write(iso_id,'(i0)') which_iso
        do which_L = 0,7
        write(L_id,'(i0)') which_L
        do which_k = 0, which_L
           do which_ph = 1,15
              write(k_id,'(i0)') which_k
              write(ph_id,'(i0)') which_ph
                if(which_k.eq.0) then
                   Ph_mult = 1                   !!!Q: apply for deformed nuclear.          
                else
                   Ph_mult = 2  
                end if
                

	  inquire(file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),exist=file_exists)    !modify18
      if(file_exists) then

      nph = nph + 1

	  open(110,file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),status='unknown',form='unformatted')       !modify18
          Read(110) WX_realomega
	      Read(110) WX_iHFB_NB                         !!!Q: read


	      if(allocated(WX_id)) deallocate(WX_id, WX_sum_nd_colm1)
	      allocate(WX_id(WX_iHFB_NB), WX_sum_nd_colm1(WX_iHFB_NB))

          if(allocated(WX_readin_select))  deallocate(WX_readin_select)
	  Allocate(WX_readin_select(WX_iHFB_NB,WX_iHFB_NB))


	  Read(110) WX_NumberofPH11                    !!!Q: read
	  Read(110) WX_id                              !!!Q: read
	  Read(110) WX_readin_select                   !!!Q: read
	  !
	  close(110)






          if(allocated(WX_NH11_E101)) deallocate(WX_NH11_E101, & 
                                                 WX_PH11_E101)
	  Allocate(WX_NH11_E101(WX_NumberofPH11),WX_PH11_E101(WX_NumberofPH11))
	  WX_NH11_E101=0
	  WX_PH11_E101=0

          !--------------------------------------------------
    !     if(allocated(WX_BNH11_E101)) deallocate(WX_BNH11_E101, & 
    !                                            WX_BPH11_E101)
	! Allocate(WX_BNH11_E101(WX_NumberofPH11),WX_BPH11_E101(WX_NumberofPH11))
	! WX_BNH11_E101=0
	! WX_BPH11_E101=0
          !--------------------------------------------------

          !----------------------------------------E101
          ! 
	  open(201,file=trim(adjustl(qq_io_path))//'QQ_binary_NH11_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),& 
                                                                status='unknown',form='unformatted')        !modify18
	  Read(201) WX_NH11_E101
	! Read(201) WX_BNH11_E101
	  close(201)
	  !
	  !
	  open(201,file=trim(adjustl(qq_io_path))//'QQ_binary_PH11_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),&
                                                                status='unknown',form='unformatted')        !modify18
	  Read(201) WX_PH11_E101
	! Read(201) WX_BPH11_E101
	  close(201)
          ! 

     WX_id_store(nph,:) = WX_id
     WX_readin_select_store(nph,:,:) = WX_readin_select
!
    do qqi = 1, size(WX_NH11_E101)
     ! WX_NH11_store(nph, qqi) = WX_NH11_E101 (qqi)
     ! WX_BNH11_store(nph,qqi) = WX_BNH11_E101(qqi)  
     ! WX_PH11_store(nph, qqi) = WX_PH11_E101(qqi)
     ! WX_BPH11_store(nph,qqi) = WX_BPH11_E101(qqi)


       WX_NH11_store1D((nph-1)*WX_NumberofPH11_maxval + qqi) = WX_NH11_E101(qqi)
     ! WX_BNH11_store1D((nph-1)*WX_NumberofPH11_maxval + qqi) = WX_BNH11_E101(qqi)
       WX_PH11_store1D((nph-1)*WX_NumberofPH11_maxval + qqi) = WX_PH11_E101(qqi)
     ! WX_BPH11_store1D((nph-1)*WX_NumberofPH11_maxval + qqi) = WX_BPH11_E101(qqi)
    enddo
!
    


      end if    ! file_exists
         end do ! which_ph
         end do ! which_k
         end do ! which_L
         end do ! which_iso


   !do qqi = 1, size(WX_NH11_E101)
   !  if(WX_NH11_store(nph, qqi).Ne. WX_NH11_E101(qqi) .and. &  
   !    WX_NH11_store(nph, qqi).Ne. WX_NH11_E101(qqi) .and.&
   !    WX_NH11_store(nph, qqi).Ne. WX_NH11_E101(qqi) .and.& 
   !    WX_NH11_store(nph, qqi).Ne. WX_NH11_E101(qqi))  then 
   !      write(50,*) 'Wrong element assignment in WX_NH11_store(nph))'
   !  endif
   !enddo
   !write(50,*) 'WX_NH11_store(nph) has same element as WX_NH11_E101'
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
        ! 
      end if			!!!Q: preparation before iter loop end.
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! end preparation
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 


      time1=0; time2=0; timei=0; timei0=0

      ! Iteration header
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,5a)') '    i','               si',&
         '               Re(S)','               Im(S)','      Time'; call writelog(st)
      write(st,'(1x,72("-"))'); call writelog(st)

      ! Main iteration loop
      do iter=0,max_iter

         ! Record the iteration (or show where we are starting from)
         write(st,'(1x, i4,a1, 2x, f15.10, 2x, f18.10, 2x, f18.10, 2x, f8.3)') &
            iter, qrpa_bbroyden, si, re_str(1), im_str(1), timei; call writelog(st)
         if (allocated(g)) then
            do ixterm=1, size(g)
               write(st,'(1x,2x,a10,2x,"Re(X)=",f16.8,2x,"Im(X)=",f16.8)') &
                   trim(g(ixterm)%label), re_str(1+ixterm), im_str(1+ixterm); call writelog(st)
            end do
         end if

         ! Exit condition
         if (si < convergence_epsilon) then
            iter_conv = iter; exit
         end if
         if (iter==max_iter) then
            iter_conv = -1; exit
         end if

         ! Start the iteration
         time1 = get_timer()
         if (iter==0) timei0 = time1

         if (abs(quench_residual_int) < 1e-10) then
            call set_val_bbm(dHqp_im,0d0)
            call set_val_bbm(dHqp_re,0d0)
         else

            ! Compute the perturbed densities in SP space
            if (bminus) then
               call triprod_bbm('n',Wp,'n',dRqp_re,'t',Wn,dRsp_re)
               call triprod_bbm('n',Wp,'n',dRqp_im,'t',Wn,dRsp_im)
            else
               call triprod_bbm('n',Wn,'n',dRqp_re,'t',Wp,dRsp_re)
               call triprod_bbm('n',Wn,'n',dRqp_im,'t',Wp,dRsp_im)
            end if

            call calc_hamiltonian(&
                                  dRsp_re%m11, dRsp_im%m11, dRsp_re%m12, dRsp_im%m12,&
                                  dRsp_re%m22, dRsp_im%m22, dRsp_re%m21, dRsp_im%m21,&
                                  dHsp_re%m11, dHsp_im%m11, dHsp_re%m12, dHsp_im%m12,&
                                  dHsp_re%m22, dHsp_im%m22, dHsp_re%m21, dHsp_im%m21 &
                                  )

            ! Transform the Hamiltonian into quasiparticle basis
            if (bminus) then
               call triprod_bbm('t',Wp,'n',dHsp_re,'n',Wn,dHqp_re)
               call triprod_bbm('t',Wp,'n',dHsp_im,'n',Wn,dHqp_im)
            else
               call triprod_bbm('t',Wn,'n',dHsp_re,'n',Wp,dHqp_re)
               call triprod_bbm('t',Wn,'n',dHsp_im,'n',Wp,dHqp_im)
            end if

            ! Quench the residual interaction
            call scalar_mult_bbm(dHqp_re,quench_residual_int)
            call scalar_mult_bbm(dHqp_im,quench_residual_int)

         end if ! End computation of the residual interaction


         !------------------------------------------QQst        modify21
         ! assign array X and Y to Ph_matX and Ph_matY to construct product 
         !----------------------------
        if(Ph_test) then

          !--------------------------------
          ! readin binfo, N11 and P11 from lpfam
          !--------------------------------
	  !
        do nph = 1, nph_rec
           if(k_store(nph) .eq. 0) then 
              Ph_mult = 1
           else
              Ph_mult = 1
           end if
     !   do which_iso=  0,1
     !   write(iso_id,'(i0)') which_iso
     !   do which_L = 0,7
     !   write(L_id,'(i0)') which_L
     !   do which_k = 0, which_L
     !      do which_ph = 1,15
     !         write(k_id,'(i0)') which_k
     !         write(ph_id,'(i0)') which_ph
     !           if(which_k.eq.0) then
     !              Ph_mult = 1                   !!!Q: apply for deformed nuclear.          
     !           else
     !              Ph_mult = 2  
     !           end if
     !           
!
	  !inquire(file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),exist=file_exists)    !modify18
     ! if(file_exists) then
     ! 
     ! nph = nph + 1
    
	 !open(110,file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),status='unknown',form='unformatted')       !modify18
     !    Read(110) WX_realomega
	 !Read(110) WX_iHFB_NB                         !!!Q: read


	  if(allocated(WX_sum_nd_colm1)) deallocate(WX_sum_nd_colm1)
	  allocate(WX_sum_nd_colm1(WX_iHFB_NB))

     !    if(allocated(WX_readin_select))  deallocate(WX_readin_select)
	 !Allocate(WX_readin_select(WX_iHFB_NB,WX_iHFB_NB))


	 !Read(110) WX_NumberofPH11                    !!!Q: read
	 !Read(110) WX_id                              !!!Q: read
	 !Read(110) WX_readin_select                   !!!Q: read
	  !
	 !close(110)

     W_omega = 1.0
     if(WX_realomega_store(nph) .lt.0.5) then
          write(50,*) iter, 'TLK: ', iso_store(nph), L_store(nph), k_store(nph), 'cycle', WX_realomega_store(nph)
          cycle
     end if
 


         !open(600,file='600_Ph_db_WX_id')
         !write(600, *) Ph_db
         !write(600, *) WX_id_store(nph,:)
         !close(600)
!-------------------------------------------------------------------
          if(allocated(WX_pnbasis_select)) deallocate(WX_pnbasis_select)
          Allocate(WX_pnbasis_select(2*WX_iHFB_NB, 2*WX_iHFB_NB))
          !
	  WX_sum_nd=0
	  WX_sum_nd_colm1=0
	  do qqj=1,WX_iHFB_NB
	      WX_sum_nd=WX_sum_nd + WX_id_store(nph, qqj)
	      WX_sum_nd_colm1(qqj)=WX_sum_nd - WX_id_store(nph, qqj)
	  end do
	  ! 
	  !
          if(allocated(WX_matN11_E101)) deallocate(WX_matN11_E101, &
                                        WX_matP11_E101)
	  Allocate(WX_matN11_E101(WX_sum_nd,WX_sum_nd),WX_matP11_E101(WX_sum_nd,WX_sum_nd))
	  WX_matN11_E101=0
	  WX_matP11_E101=0
          !
          ! 
          !-----------------------------------------------------------
          if(allocated(WX_matBN11_E101)) deallocate(WX_matBN11_E101, &
                                        WX_matBP11_E101)
	  Allocate(WX_matBN11_E101(WX_sum_nd,WX_sum_nd),WX_matBP11_E101(WX_sum_nd,WX_sum_nd))
	  WX_matBN11_E101=0
	  WX_matBP11_E101=0
          !-----------------------------------------------------------
          !
          !
          if(allocated(WX_matpnbasisN11_E101)) deallocate(WX_matpnbasisN11_E101,&
                                                WX_matpnbasisP11_E101) 
	  Allocate(WX_matpnbasisN11_E101(2*WX_sum_nd,2*WX_sum_nd),WX_matpnbasisP11_E101(2*WX_sum_nd,2*WX_sum_nd))
	  WX_matpnbasisN11_E101=0
	  WX_matpnbasisP11_E101=0


!---------------------------------------------------------------------------------



          write(50,*) iter, trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)), & 
                 Ph_sum_db, 2*WX_sum_nd, Ph_mult, WX_realomega_store(nph) 

     !    if(allocated(WX_NH11_E101)) deallocate(WX_NH11_E101, & 
     !                                           WX_PH11_E101)
	 !Allocate(WX_NH11_E101(WX_NumberofPH11),WX_PH11_E101(WX_NumberofPH11))
	 !WX_NH11_E101=0
	 !WX_PH11_E101=0

     !    !--------------------------------------------------
     !    if(allocated(WX_BNH11_E101)) deallocate(WX_BNH11_E101, & 
     !                                           WX_BPH11_E101)
	 !Allocate(WX_BNH11_E101(WX_NumberofPH11),WX_BPH11_E101(WX_NumberofPH11))
	 !WX_BNH11_E101=0
	 !WX_BPH11_E101=0
     !    !--------------------------------------------------


         !Call WgreenX_func(Ph_sum_db, Ep, En, omega, WX_realomega_store(nph), width, WgreenX) 
          Call WgreenX_func(Ph_sum_db, Ep, En, real_eqrpa, WX_realomega_store(nph), imag_eqrpa, WgreenX) 
         !Call WgreenY_func(Ph_sum_db, Ep, En, omega, WX_realomega_store(nph), width, WgreenY) 
          Call WgreenY_func(Ph_sum_db, Ep, En, real_eqrpa, WX_realomega_store(nph), imag_eqrpa, WgreenY) 
 
          do qqi = 1, 2*WX_iHFB_NB
             do qqj = 1, 2*WX_iHFB_NB
                if ((qqi.le.WX_iHFB_NB).and.(qqj.le.WX_iHFB_NB)) then
                        WX_pnbasis_select(qqi, qqj) = WX_readin_select_store(nph, qqi, qqj)
                end if
                !
                if ((qqi.le.WX_iHFB_NB).and.(qqj.gt.WX_iHFB_NB)) then
                        WX_pnbasis_select(qqi, qqj) = WX_readin_select_store(nph, qqi, qqj-WX_iHFB_NB)
                end if
                !
                if ((qqi.gt.WX_iHFB_NB).and.(qqj.le.WX_iHFB_NB)) then
                        WX_pnbasis_select(qqi, qqj) = WX_readin_select_store(nph, qqi-WX_iHFB_NB, qqj)
                end if
                !
                if ((qqi.gt.WX_iHFB_NB).and.(qqj.gt.WX_iHFB_NB)) then
                        WX_pnbasis_select(qqi, qqj) = WX_readin_select_store(nph, qqi-WX_iHFB_NB, qqj-WX_iHFB_NB)
                end if
             end do
          end do
 


          !----------------------------------------E101
          ! 
	  !open(201,file=trim(adjustl(qq_io_path))//'QQ_binary_NH11_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),& 
     !                                                           status='unknown',form='unformatted')        !modify18
	  !Read(201) WX_NH11_E101
	  !Read(201) WX_BNH11_E101
	  !close(201)
	  !!
	  !!
	  !open(201,file=trim(adjustl(qq_io_path))//'QQ_binary_PH11_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),&
     !                                                           status='unknown',form='unformatted')        !modify18
	  !Read(201) WX_PH11_E101
	  !Read(201) WX_BPH11_E101
	  !close(201)
          ! 
	  !-------------------------
	  ! asssign NH11 to matN11
	  ! asssign NP11 to matP11
	  !-------------------------
	  !
	  ct1=1
	  ct3=0
	  inbcc=0
	  inbrc=0
	  inbrr=0
	  do qqj=1,WX_iHFB_NB
	     do qqi=1,WX_iHFB_NB
		if(WX_readin_select_store(nph, qqi,qqj)==1)  then
		   inbcc=WX_sum_nd_colm1(qqj)+1
		   inbrc=WX_sum_nd_colm1(qqi)
		   ct3=ct3+WX_id_store(nph, qqj)*WX_id_store(nph, qqi)
		   do ct2=ct1, ct3 
		      inbrc=inbrc+1
		      inbrr=inbrr+1
                      !----------------------------E101
		      WX_matN11_E101(inbrc, inbcc)=WX_NH11_store1D((nph-1)*WX_NumberofPH11_maxval+ct2)              
		     !WX_matBN11_E101(inbrc, inbcc)=WX_BNH11_store1D((nph-1)*WX_NumberofPH11_maxval+ct2)              
                      !
		      WX_matP11_E101(inbrc, inbcc)=WX_PH11_store1D((nph-1)*WX_NumberofPH11_maxval+ct2)              
		     !WX_matBP11_E101(inbrc, inbcc)=WX_BPH11_store1D((nph-1)*WX_NumberofPH11_maxval+ct2)              
                      !
		      if(inbrr.eq. WX_id_store(nph, qqi)) then
			 inbcc=inbcc+1        !!!Q: col increase by 1
			 inbrc=WX_sum_nd_colm1(qqi)   !!!Q: reset inbrc   
			 inbrr=0              !!!Q: reset inbrr
		      end if
		   end do
		   ct1=ct1+WX_id_store(nph, qqj)*WX_id_store(nph, qqi)
		end if
	     end do
	  end do
          !
          !----------------------------------------------------------
	

          !-----------------------------------QQst      modify22
          do qqj = 1, Ph_sum_db   
             do qqi = 1, Ph_sum_db
                ! (row positive Omega, col positive Omega)
                if((qqi.Le.WX_sum_nd) .and.(qqj.Le.WX_sum_nd)) then 
                        !-----------------------------E101
                        WX_matpnbasisN11_E101(qqi,qqj) = &
                             !0.5*(WX_matN11_E101(qqi,qqj)+WX_matBN11_E101(qqi,qqj))
                              real(WX_matN11_E101(qqi,qqj))*2
                        WX_matpnbasisP11_E101(qqi,qqj) = &
                             !0.5*(WX_matP11_E101(qqi,qqj)+WX_matBP11_E101(qqi,qqj))
                              real(WX_matP11_E101(qqi,qqj))*2
                                 
                end if
                ! (row positive Omega, col negative Omega) 
                if((qqi.Le.WX_sum_nd) .and.(qqj.gt.WX_sum_nd)) then 
                        !-----------------------------E101
                        WX_matpnbasisN11_E101(qqi,qqj) = &
                            ! 0.5*cunit*(WX_matN11_E101(qqi,qqj-WX_sum_nd)-WX_matBN11_E101(qqi,qqj-WX_sum_nd))
                             +aimag(WX_matN11_E101(qqi,qqj-WX_sum_nd))*2
                        WX_matpnbasisP11_E101(qqi,qqj) = &
                            ! 0.5*cunit*(WX_matP11_E101(qqi,qqj-WX_sum_nd)-WX_matBP11_E101(qqi,qqj-WX_sum_nd))
                             +aimag(WX_matP11_E101(qqi,qqj-WX_sum_nd))*2
                end if
                ! (row negative Omega, col positive Omega) 
                if((qqi.gt.WX_sum_nd) .and.(qqj.Le.WX_sum_nd)) then 
                        !-----------------------------E101
                        WX_matpnbasisN11_E101(qqi,qqj) = &
                            !-0.5*cunit*(WX_matN11_E101(qqi-WX_sum_nd,qqj)-WX_matBN11_E101(qqi-WX_sum_nd,qqj))
                             -aimag(WX_matN11_E101(qqi-WX_sum_nd,qqj))*2
                        WX_matpnbasisP11_E101(qqi,qqj) = &
                            !-0.5*cunit*(WX_matP11_E101(qqi-WX_sum_nd,qqj)-WX_matBP11_E101(qqi-WX_sum_nd,qqj))
                             -aimag(WX_matP11_E101(qqi-WX_sum_nd,qqj))*2
                end if
                ! (row negative Omega, col negative Omega) 
                if((qqi.gt.WX_sum_nd) .and.(qqj.gt.WX_sum_nd)) then 
                        !-----------------------------E101
                        WX_matpnbasisN11_E101(qqi,qqj) = &
                            ! 0.5*(WX_matN11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd)+WX_matBN11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd))
                              real(WX_matN11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd))*2
                        WX_matpnbasisP11_E101(qqi,qqj) = &
                            ! 0.5*(WX_matP11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd)+WX_matBP11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd))
                              real(WX_matP11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd))*2
                end if
             end do
          end do
          !----------------------------------QQend             


         !=================================================
         !=================================================
         !=================================================

         !Call WX_cal          
          call WX_cal(nxy,Ph_nb, Ph_db, Ph_readin_select, &
                     Ph_sum_db,Ph_sum_db_colm1, Ph_NumberofX, WX_sum_nd, & 
                    !ReX%elem, ImX%elem, ReY%elem, ImY%elem, &
                     dHqp_re%m12%elem, dHqp_im%m12%elem,dHqp_re%m21%elem, dHqp_im%m21%elem,  &  
                     WX_matpnbasisN11_E101, WX_matpnbasisP11_E101, &
                     Ph_ReWX, Ph_ImWX, Ph_ReWY, Ph_ImWY, &
                     WgreenX, WgreenY,&
                     WX_iHFB_NB,WX_pnbasis_select,X_select,Y_select)       

            !h20re%elem = h20re%elem - Ph_mult*Ph_ReWX*W_omega 
             dHqp_re%m12%elem = dHqp_re%m12%elem - Ph_mult*Ph_ReWX*W_omega
            !h20Im%elem = h20Im%elem - Ph_mult*Ph_ImWX*W_omega 
             dHqp_im%m12%elem = dHqp_im%m12%elem - Ph_mult*Ph_ReWX*W_omega
             !
            !h02re%elem = h02re%elem - Ph_mult*Ph_ReWY*W_omega 
             dHqp_re%m21%elem = dHqp_re%m21%elem + Ph_mult*Ph_ReWX*W_omega
            !h02Im%elem = h02Im%elem + Ph_mult*Ph_ImWY*W_omega    !!!Q: need '-' as (WY*)* conjg.
             dHqp_im%m21%elem = dHqp_im%m21%elem - Ph_mult*Ph_ReWX*W_omega

        !write(*,*) nxy, sizeof(Fqp), sizeof(dHqp_re), sizeof(dRqp_re)
        !write(*,*) nxy, size(dHqp_re%m11%elem),size(dHqp_re%m12%elem),size(dHqp_re%m21%elem),size(dHqp_re%m22%elem)
        !call effective_WX(nxy, dRqp_re%m12%elem, dRqp_im%m12%elem, dRqp_re%m21%elem,dRqp_im%m21%elem)





           ! h20re%elem = h20re%elem + 5.000*Ph_ReWX 
           ! h20Im%elem = h20Im%elem + 5.000*Ph_ImWX 
           ! !
           ! h02re%elem = h02re%elem + 5.000*Ph_ReWY 
           ! h02Im%elem = h02Im%elem - 5.000*Ph_ImWY    !!!Q: need '-' as (WY*)* conjg.
 
 

      !end if    ! file_exists
      !   end do ! which_ph
      !   end do ! which_k
      !   end do ! which_L
      !   end do ! which_iso
      end do    ! nph
        end if  ! Ph_test
      
        !------------------------------------------QQend 




         ! Add the external field to perturbed hamiltonian (assuming real Fqp)
         call add_bbm(Fqp, dHqp_re)

         ! Solve the new X/Y, P/Q
         call complex_emult_bbm(Greens_re, Greens_im, dHqp_re, dHqp_im,&
                                         dRqp_re, dRqp_im)

        !write(*,*) nxy, sizeof(Fqp), sizeof(dHqp_re), sizeof(dRqp_re)
         write(*,*) nxy, size(dHqp_re%m11%elem),size(dHqp_re%m12%elem),size(dHqp_re%m21%elem),size(dHqp_re%m22%elem)
         call effective_WX(nxy, dRqp_re%m12%elem, dRqp_im%m12%elem, dRqp_re%m21%elem,dRqp_im%m21%elem)


         if (use_diagonal_blocks) then
            call emult_bbm(T, dRqp_re, dRqp_re)
            call emult_bbm(T, dRqp_im, dRqp_im)
         end if

         ! Use Broyden mixing to stabilize convergence
         ! T > 0 / odd-A have twice the number of vectors to handle P and Q
         call qrpa_broyden(iter, dRqp_re, dRqp_im)

         ! Strength function
         call contract_bbm(Fqp,dRqp_re,re_s)
         call contract_bbm(Fqp,dRqp_im,im_s)
         re_str(1) = -re_s/pi
         im_str(1) = -im_s/pi

         ! Cross-terms
         if (allocated(g)) then
            do ixterm=1, size(g)
               call contract_bbm(Gqp(ixterm),dRqp_re,re_s)
               call contract_bbm(Gqp(ixterm),dRqp_im,im_s)
               re_str(1+ixterm) = -re_s/pi
               im_str(1+ixterm) = -im_s/pi
            end do
         end if

         ! Record the iteration time
         time2 = get_timer()
         timei = (time2-time1)

      end do

      write(st,'(1x,72("-"))'); call writelog(st)
      if (iter_conv >= 0 ) then
         write(st,'("  *  FAM iteration converged after   ",i4," steps.  si = ",ES12.5)') iter,si; call writelog(st)
      else
         write(st,'("  *  FAM iteration interrupted after ",i4," steps.  si = ",ES12.5)') iter,si; call writelog(st)
      end if
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,a,es12.3,a)') 'Total CPU time = ',(time2-timei0)/60.0_dp,' minutes' ; call writelog(st)
      call writelog("")

   end subroutine ifam

   !---------------------------------------------------------------------------
   ! Allocates and initializes the data structures used by pnfam_solve.
   !---------------------------------------------------------------------------
   subroutine init_pnfam_solver
      use pnfam_constants, only : IT_NEUTRON, IT_PROTON
      use pnfam_broyden, only : qrpa_broin, qrpa_broout, nbroy, broyden_history_size

      implicit none

      integer :: ixterm, i
      logical :: a11,a12,a21,a22

      !---------------------------------------------------------------------------
      ! General initializations
      !---------------------------------------------------------------------------
      ! The complex strength function (array [strength, xterms])
      if (allocated(re_str)) deallocate(re_str, im_str)
      allocate(re_str(1:nxterms+1), im_str(1:nxterms+1))
      re_str(:) = 0; im_str(:) = 0

      ! While the various FAM matrices may have different non-trivial block structures
      ! they still share the same number of elements in the non-trivial blocks
      nxy = size(f%mat%elem)
      bminus = f%beta_minus

      if (ft_active .or. hfb_blo_active) then
         use_diagonal_blocks = .true.
         nbroy = nxy*8
      else
         use_diagonal_blocks = .false.
         nbroy = nxy*4
      end if

      ! Broyden arrays
      if (allocated(qrpa_broin)) deallocate(qrpa_broin,qrpa_broout)
      allocate(qrpa_broin(nbroy),qrpa_broout(nbroy))
      qrpa_broin = 0; qrpa_broout=0
      si = 1
      ! Don't do broyden if no residual interaction, should just need 1 iteration
      if (abs(quench_residual_int) < 1e-10) then
         broyden_history_size = -1
      end if

      !---------------------------------------------------------------------------
      ! Allocate and initialize big-blockmatrix structure
      !---------------------------------------------------------------------------
      ! Which blocks are active
      a11 = .false.; a12 = .true.
      a21 = .true. ; a22 = .false.
      if (use_diagonal_blocks) then
         a11 = .true.; a22 = .true.
      end if

      ! Single particle external field (F_ij*(a^\dagger a) --> m11)
      ! Represents.. : [ f_pn, 0 ; 0, 0 ], for 2bc w/ pairing: [f11/2, f12; -f12*/2, -f11*]
      ! Stores...... : [ f_pn, _ ; _, _ ], for 2bc w/ pairing: [f11/2, f12;  f12/2,   f11 ]
      Fsp%m11 = f%mat
      ! THIS DOES NOT WORK B/C IT SPOILS THE 1 BLOCK PER ROW+COL STRUCTURE IN FQP
      if (allocated(f%mat12%elem)) then
         Fsp%m12 = f%mat12
         Fsp%m21 = f%mat12
         Fsp%m22 = f%mat
         Fsp%m11%elem = f%mat%elem*0.5_dp
         Fsp%m22%elem = f%mat%elem*0.5_dp
         call set_sign_bbm(Fsp, 1d0,1d0,-1d0,-1d0)
      end if
      if (allocated(g)) then
         if (allocated(Gsp)) deallocate(Gsp)
         allocate(Gsp(size(g)), Gqp(size(g)))
         do ixterm=1, size(g)
            Gsp(ixterm)%m11 = g(ixterm)%mat
         end do
      end if

      ! Bogoliubov transform (real matrices)
      ! Represents.. : [ U, V* ; V, U* ]
      ! Stores...... : [ U, V  ; V, U  ]
      Wn%m11 = Un; Wn%m12 = Vn
      Wn%m21 = Vn; Wn%m22 = Un
      Wp%m11 = Up; Wp%m12 = Vp
      Wp%m21 = Vp; Wp%m22 = Up

      ! External fields in QP basis
      ! Represents.. : [ F11, F20 ; -F02, -F11t ]
      ! Stores...... : [ F11, F20 ;  F02,  F11t ]
      call allocate_bbm(Fqp, a11,a12,a21,a22, nxy)
      call set_sign_bbm(Fqp, 1d0,1d0,-1d0,-1d0)
      if (allocated(g)) then
         do ixterm=1, size(g)
            call allocate_bbm(Gqp(ixterm), a11,a12,a21,a22, nxy)
            call set_sign_bbm(Gqp(ixterm), 1d0,1d0,-1d0,-1d0)
         end do
      end if

      ! Perturbed hamiltonian in QP basis
      ! Represents.. : [ H11, H20 ; -H02, -H11t ]
      ! Stores...... : [ H11, H20 ;  H02,  H11t ]
      call allocate_bbm(dHqp_re, a11,a12,a21,a22, nxy)
      call allocate_bbm(dHqp_im, a11,a12,a21,a22, nxy); dHqp_im%imag=.true.
      call set_sign_bbm(dHqp_re, 1d0,1d0,-1d0,-1d0)
      call set_sign_bbm(dHqp_im, 1d0,1d0,-1d0,-1d0)

      ! Perturbed generalized density in QP basis
      ! Represents.. : [ P, X ; -Y, -Q ]
      ! Stores...... : [ P, X ;  Y,  Q ]
      call allocate_bbm(dRqp_re, a11,a12,a21,a22, nxy)
      call allocate_bbm(dRqp_im, a11,a12,a21,a22, nxy); dRqp_im%imag=.true.
      call set_sign_bbm(dRqp_re, 1d0,1d0,-1d0,-1d0)
      call set_sign_bbm(dRqp_im, 1d0,1d0,-1d0,-1d0)

      ! FAM greens function: g = -1/(E1 +/- E2 +/- (re(eqrpa)+im(eqrpa)*i))
      ! Represents.. : diag( gP, gX, gY, gQ )
      ! Stores...... : [ gP, gX ; gY, gQ ]
      call allocate_bbm(Greens_re, a11,a12,a21,a22, nxy)
      call allocate_bbm(Greens_im, a11,a12,a21,a22, nxy); Greens_im%imag=.true.

      ! Statistical occupation matrix: 1ff=(1-f1-f2) and ff=(f2-f1)
      ! Represents.. : diag( T1ff Tff Tff T1ff )
      ! Stores...... : [ T1ff, Tff ; Tff, T1ff ]
      if (use_diagonal_blocks) call allocate_bbm(T, a11,a12,a21,a22, nxy)

      ! Perturbed hamiltonian in SP basis
      ! Represents.. : [ h_pn, Delta(+)_pn ; -Delta(-)*_pn, -(h^T)_pn ]
      ! Stores...... : [ h_pn, Delta(+)_pn ;  Delta(-) _pn,   h_np    ]
      call allocate_bbm(dHsp_re, .true.,.true.,.true.,.true., nxy)
      call allocate_bbm(dHsp_im, .true.,.true.,.true.,.true., nxy); dHsp_im%imag=.true.
      call  set_sign_bbm(dHsp_re,+1d0,+1d0,-1d0,-1d0)
      call set_trans_bbm(dHsp_re,'n', 'n', 'n', 't')
      call  set_sign_bbm(dHsp_im,+1d0,+1d0,+1d0,-1d0)
      call set_trans_bbm(dHsp_im,'n', 'n', 'n', 't')

      ! Perturbed density in SP basis (Note kap is antisymmetric so np<-->pn gives a minus sign)
      ! Represents.. : [ rho_pn,  kap(+)_pn ; -kap(-)*_pn, -(rho^T)_pn ]
      ! Stores...... : [ rho_pn,  kap(+)_np ;  kap(-) _np,   rho_np    ]
      !            ?=? [ rho_pn, rhotilde(+)_pn ; rhotilde(-)_pn, rho_np ]
      call allocate_bbm(dRsp_re, .true.,.true.,.true.,.true., nxy)
      call allocate_bbm(dRsp_im, .true.,.true.,.true.,.true., nxy); dRsp_im%imag=.true.
      call  set_sign_bbm(dRsp_re,+1d0,-1d0,+1d0,-1d0)
      call set_trans_bbm(dRsp_re,'n', 't', 't', 't')
      call  set_sign_bbm(dRsp_im,+1d0,-1d0,-1d0,-1d0)
      call set_trans_bbm(dRsp_im,'n', 't', 't', 't')

      !---------------------------------------------------------------------------
      ! Transform external field to QP basis
      !---------------------------------------------------------------------------
      if (bminus) then
         call triprod_bbm('t',Wp,'n',Fsp,'n',Wn,Fqp)
      else
         call triprod_bbm('t',Wn,'n',Fsp,'n',Wp,Fqp)
      end if
      if (allocated(g)) then
         do ixterm=1, size(g)
            if (bminus) then
               call triprod_bbm('t',Wp,'n',Gsp(ixterm),'n',Wn,Gqp(ixterm))
            else
               call triprod_bbm('t',Wn,'n',Gsp(ixterm),'n',Wp,Gqp(ixterm))
            end if
         end do
      end if

      !---------------------------------------------------------------------------
      ! Initialize individual block structures based on external field
      ! (For arrays whose block structure is not determined by routine triprod)
      !---------------------------------------------------------------------------
      ! dRqp has same structure as Fqp
      call copy_block_structure_bbm(Fqp, dRqp_re)
      call copy_block_structure_bbm(Fqp, dRqp_im)

      ! dHqp has same structure as Fqp
      call copy_block_structure_bbm(Fqp, dHqp_re)
      call copy_block_structure_bbm(Fqp, dHqp_im)

      ! Greens has same structure as dRqp
      call copy_block_structure_bbm(dRqp_re, Greens_re)
      call copy_block_structure_bbm(dRqp_im, Greens_im)

      ! Statmatrix has same structure as Fqp
      if (use_diagonal_blocks) call copy_block_structure_bbm(Fqp, T)

      ! dh_pn has same structure as Fsp
      call copy_block_structure(Fsp%m11, dHsp_re%m11)
      call copy_block_structure(Fsp%m11, dHsp_im%m11)
      ! dDelta(+)_pn has same structure as X
      call copy_block_structure(dRqp_re%m12, dHsp_re%m12)
      call copy_block_structure(dRqp_im%m12, dHsp_im%m12)
      ! dDelta(-)_pn has same structure as Y
      call copy_block_structure(dRqp_re%m21, dHsp_re%m21)
      call copy_block_structure(dRqp_im%m21, dHsp_im%m21)
      ! dh_np has same structure as Fsp
      call copy_block_structure(Fsp%m11, dHsp_re%m22)
      call copy_block_structure(Fsp%m11, dHsp_im%m22)

      !---------------------------------------------------------------------------
      ! Calculate 2QP quantities
      !---------------------------------------------------------------------------
      ! Greens function: -1/(E1 +/- E2 +/- (re(eqrpa)+im(eqrpa)*i))
      ! Calculate the 2QP matrix elements M_{ij} = a*( b*f1_i + c*f2_j + e )**d
      if (bminus) then
         call matrix_2qp(-1,+1,+1,-1,-real_eqrpa,-imag_eqrpa,.true., Ep,En,Greens_re%m12)
         call matrix_2qp(-1,+1,+1,-1,-real_eqrpa,-imag_eqrpa,.false.,Ep,En,Greens_im%m12)
         call matrix_2qp(-1,+1,+1,-1,+real_eqrpa,+imag_eqrpa,.true., Ep,En,Greens_re%m21)
         call matrix_2qp(-1,+1,+1,-1,+real_eqrpa,+imag_eqrpa,.false.,Ep,En,Greens_im%m21)
         if (use_diagonal_blocks) then
         call matrix_2qp(-1,+1,-1,-1,-real_eqrpa,-imag_eqrpa,.true., Ep,En,Greens_re%m11)
         call matrix_2qp(-1,+1,-1,-1,-real_eqrpa,-imag_eqrpa,.false.,Ep,En,Greens_im%m11)
         call matrix_2qp(-1,+1,-1,-1,+real_eqrpa,+imag_eqrpa,.true., Ep,En,Greens_re%m22)
         call matrix_2qp(-1,+1,-1,-1,+real_eqrpa,+imag_eqrpa,.false.,Ep,En,Greens_im%m22)
         end if
      else
         call matrix_2qp(-1,+1,+1,-1,-real_eqrpa,-imag_eqrpa,.true., En,Ep,Greens_re%m12)
         call matrix_2qp(-1,+1,+1,-1,-real_eqrpa,-imag_eqrpa,.false.,En,Ep,Greens_im%m12)
         call matrix_2qp(-1,+1,+1,-1,+real_eqrpa,+imag_eqrpa,.true., En,Ep,Greens_re%m21)
         call matrix_2qp(-1,+1,+1,-1,+real_eqrpa,+imag_eqrpa,.false.,En,Ep,Greens_im%m21)
         if (use_diagonal_blocks) then
         call matrix_2qp(-1,+1,-1,-1,-real_eqrpa,-imag_eqrpa,.true., En,Ep,Greens_re%m11)
         call matrix_2qp(-1,+1,-1,-1,-real_eqrpa,-imag_eqrpa,.false.,En,Ep,Greens_im%m11)
         call matrix_2qp(-1,+1,-1,-1,+real_eqrpa,+imag_eqrpa,.true., En,Ep,Greens_re%m22)
         call matrix_2qp(-1,+1,-1,-1,+real_eqrpa,+imag_eqrpa,.false.,En,Ep,Greens_im%m22)
         end if
      end if

      ! Statistical occupation matrix: (1-f1-f2) and (f2-f1)
      if (use_diagonal_blocks) then
         call set_val_bbm(T,0d0)
         if (bminus) then
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m12)
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m21)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m11)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m22)
         else
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m12)
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m21)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m11)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m22)
         end if
      end if

   end subroutine init_pnfam_solver


   subroutine deallocate_pnfam_solver
      implicit none
      integer :: ixterm

      ! External fields in QP basis
      call deallocate_bbm(Fqp)
      if (allocated(g)) then
         do ixterm=1, size(g)
            call deallocate_bbm(Gqp(ixterm))
         end do
      end if

      ! Perturbed hamiltonian in QP basis
      call deallocate_bbm(dHqp_re)
      call deallocate_bbm(dHqp_im)

      ! Perturbed generalized density in QP basis
      call deallocate_bbm(dRqp_re)
      call deallocate_bbm(dRqp_im)

      ! FAM greens function: g = -1/(E1 +/- E2 +/- (re(eqrpa)+im(eqrpa)*i))
      call deallocate_bbm(Greens_re)
      call deallocate_bbm(Greens_im)

      ! Statistical occupation matrix: 1ff=(1-f1-f2) and ff=(f2-f1)
      if (use_diagonal_blocks) call deallocate_bbm(T)

      ! Perturbed hamiltonian in SP basis
      call deallocate_bbm(dHsp_re)
      call deallocate_bbm(dHsp_im)

      ! Perturbed density in SP basis (Note kap is antisymmetric so np<-->pn gives a minus sign)
      call deallocate_bbm(dRsp_re)
      call deallocate_bbm(dRsp_im)

   end subroutine deallocate_pnfam_solver

   !----------------------------------------------------------------------------
   ! Calculate the 2QP matrix elements M_{ij} = a*( b*f1_i + c*f2_j + e )**d
   ! Input arrays are assumed size(f) = dqp (e.g. Ep, En, qp_fp, qp_fn)
   ! - Greens(X): a= -1, b= +1, c= +1, d= -1, e= -(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - Greens(Y): a= -1, b= +1, c= +1, d= -1, e= +(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - Greens(P): a= -1, b= +1, c= -1, d= -1, e= -(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - Greens(Q): a= -1, b= +1, c= -1, d= -1, e= +(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - T(X/Y)   : a= +1, b= -1, c= -1, d= +1, e= +1,        f1= qp_fp/qp_fn, f2= qp_fn/qp_fp
   ! - T(P/Q)   : a= +1, b= +1, c= -1, d= +1, e=  0,        f1= qp_fp/qp_fn, f2= qp_fn/qp_fp
   !----------------------------------------------------------------------------
   subroutine matrix_2qp(a,b,c,d,ree,ime,re,f1,f2,mat)

      implicit none

      integer, intent(in) :: a,b,c,d
      real(dp), intent(in) :: f1(:), f2(:), ree, ime
      logical, intent(in) :: re
      type(blockmatrix), intent(inout) :: mat

      logical  :: l1, l2
      integer  :: ibr, ibc, ipt, i1, i2, i1a, i1b, i2a, i2b
      complex(dp) :: aux

      mat%elem = 0

      ipt = 1
      do ibr=1, nb
         ibc = mat%ir2c(ibr)
         if (ibc == 0) cycle
            i1a = isstart(ibr) ; i1b = i1a + db(ibr) - 1
            i2a = isstart(ibc) ; i2b = i2a + db(ibc) - 1
            do i2 = i2a,i2b ! true column index
            do i1 = i1a,i1b ! true row index
               aux = a*(b*f1(i1) + c*f2(i2) + cmplx(ree,ime,dp))**d
               if (re) then
                  mat%elem(ipt) = real(aux,dp)
               else
                  mat%elem(ipt) = aimag(aux)
               end if
               ipt = ipt + 1
            end do
            end do
      end do

   end subroutine matrix_2qp

   ! Get the external field objects
   ! (depends on namelist and basis wfs)
   !
   ! init_external_field:
   ! sets: label, beta_minus, parity, K, rank(J)
   ! calls: init_fam_mapping(depends on blockmatrix, K, parity, Kpi-of-blocks)
   ! calls: ext_field_operator (calcs operator, depends on basis, block structure)
   !
   ! setup_cross_terms just calls init_ext_field for the various cross terms
   !----------------------------------------------
   subroutine setup_extfield()
      use pnfam_constants, only : get_timer
      use pnfam_extfield, only : init_external_field, setup_crossterms, set_use_2bc
      use pnfam_storage, only : fam_io
      use pnfam_extfield_2bc, only : effective_2bc_extfield, debug
      use type_extfield_2bc, only : init_extfield_2bc_type
      use type_gamdel_2bc, only : calc_totgamdel
      implicit none
      integer :: iexit, i
      real(dp), allocatable, dimension(:) :: elem_tmp
      real(dp) :: time1, time2

      iexit=1

      ! Calculate 1 body external field
      call init_external_field(beta_type=beta_type, label=operator_name, k=operator_k, op=f)
      nxterms = 0
      if (compute_crossterms) then
         call setup_crossterms(op=f, crossterms=g, n=nxterms)
      end if

      ! For production runs (no debug) we only calculate/read the 2BC Yukawa part.
      ! The 1BC and 2BC contact part are calculated using the 1BC code and added in.
      ! use_2bc code is:
      ! * 1st digit: 1=1BC+2BC, 2=2BC only
      ! * 2nd digit: 1=Full FAM, 2=SNM+LDA, 3=ASNM+LDA, 4=DME(exc,dir=0), 5=DME(exc)+FAM(dir)
      ! * 3rd digit: 1=Gamma, 2=Delta, 3=Gamma+Delta
      if (two_body_current_mode /= 0) then
         ! Set the 2BC part of the external field type
         call set_use_2bc(f, two_body_current_mode)

         ! Store the 1-body part (if doing 1BC+2BC)
         if (allocated(elem_tmp)) deallocate(elem_tmp)
         allocate(elem_tmp(size(f%mat%elem)))
         elem_tmp = 0
         if (f%use_2bc(1) == 1) then
            elem_tmp = -f%mat%elem ! Park Sign Convention: -sigma tau + 2BC
         end if

         ! Calculate 2-body contact part using (sigma tau rho)  (if not debug).
         ! Also calculates LDA and DME contributions (which have the same contact part as FAM)
         if (abs(debug)<=1) then
            call init_external_field(beta_type=beta_type, label=operator_name, k=operator_k, op=f,&
                    use_2bc=two_body_current_mode)
            do i=1,size(f%mat%elem)
               elem_tmp(i) = elem_tmp(i) + f%mat%elem(i)
            end do
 
            ! If mode needs no FAM calc, we are done: F = -sigma tau + sigma tau f(rho)
            if (f%use_2bc(2) /= 1 .and. f%use_2bc(2) /= 5) then
               f%mat%elem = elem_tmp
               return
            end if
         end if

         ! Finally, read in the Yukawa part or calculate and store it
         f%mat%elem = 0
         call fam_io(-1,iexit)
         if (iexit > 0) then 
            write(st,'(A)') "  * Calculating 2BC matrix elements..."
            call writelog(st)

            time1 = get_timer()
            call effective_2bc_extfield(f)
            time2 = get_timer()

            write(st,'(A, ES12.3)') "  * Two-body current wall time (min): ", (time2-time1)/60.0_dp
            call writelog(st)
            call writelog("")

            ! Use only the direct part of FAM with DME
            if (f%use_2bc(2) == 5) then
                call calc_totgamdel(f%gam, f%mat%elem, 'd')
            end if

            call fam_io(-2,iexit) ! deallocates f%gam/del
         end if

         ! Add 2-body + 1-body and/or contact part
         do i=1,size(f%mat%elem)
            f%mat%elem(i) = f%mat%elem(i) + elem_tmp(i)
         end do

         ! Pairing only - NB I have not checked the contact parts contribution to pairing is correct
         if (f%use_2bc(3) == 2) then
            f%mat12%elem = f%mat%elem
            call deallocate_blockmatrix(f%mat)
         end if

      end if

   end subroutine setup_extfield


   ! NB: Fsp, dRqp, Fqp, imag_eqrpa are all globals in this module
   subroutine transition_density
      use hfb_solution
      use pnfam_constants, only : pi, IT_NEUTRON, IT_PROTON
      use type_blockmatrix
      use type_extfield_2bc, only : init_extfield_2bc_type, caux, cd
      use pnfam_extfield, only : tbc_nmlda_da1, dme_exc
      implicit none
      type(bigblockmatrix) :: dRsp_QRPA, dRqp_QRPA, dR_tmp
      real(dp) :: prefactor, N, Rij, Fij, Sij, amplitude, integral_td1, integral_td2
 
      integer :: ipt, ibx1, ibx2, nd1, nd2, i2, ix2, i1, ix1
      integer :: xl1, xs1, xl2, xs2

      real(dp) :: td(nghl), wf_1(nghl), wf_2(nghl)
      real(dp) :: rhop(nghl), rhon(nghl), rho_fac(nghl)
      character(len=1200) :: denfile, stde

      ! Redundant, but init the couplings
      call init_extfield_2bc_type(.false.)

      ! Transition ampltiude from QP basis solution:
      !   * -1/pi Im(S) = 1/(hw*pi) |<i|F|0>|^2
      !   * |<i|F|0>|^2 = -hw Im(S)
      call contract_bbm(Fqp, dRqp_im, amplitude)
      amplitude = -imag_eqrpa*amplitude ! +/- |<i| F |0>|^2
      amplitude = sqrt(abs(amplitude))  !     |<i| F |0>|

      ! Assume QRPA matrices are real, then:
      !   * Im(X_FAM) = -1/hw (X_QRPA) <i|F|0>
      !   * X_QRPA = (-hw/<i|F|0>) Im(X_FAM)
      dRqp_QRPA = dRqp_im    ! Make a copy
      dRqp_QRPA%imag=.false. ! Declare real
      prefactor = -imag_eqrpa/amplitude
      call scalar_mult_bbm(dRqp_QRPA, prefactor)

      ! Check normalization to see if we are at an eigenvalue:
      !   * N = \sum_{uv} |X_{uv}|^2 - |Y_{uv}|^2
      dR_tmp = dRqp_QRPA
      if (allocated(dR_tmp%m21%elem)) dR_tmp%m21%elem = -dR_tmp%m21%elem
      if (allocated(dR_tmp%m22%elem)) dR_tmp%m22%elem = -dR_tmp%m22%elem
      call contract_bbm(dRqp_QRPA, dR_tmp, N)

      ! Convert X/Y_QRPA to single particle basis
      dRsp_QRPA = dRsp_im     ! Make a copy
      dRsp_QRPA%imag = .false.! Declare real
      if (bminus) then
         call triprod_bbm('n',Wp,'n',dRqp_QRPA,'t',Wn, dRsp_QRPA)
      else
         call triprod_bbm('n',Wn,'n',dRqp_QRPA,'t',Wp, dRsp_QRPA)
      end if

      ! Calculate density of the transition amplitude
      !   * Since Fsp = (f11, 0; 0, 0), only need to contract with drho_QRPA
      !   * TD = \sum_{ij} rho_{ij} f11_{ij} phi_{i} phi_{j} ??? Also what about spin and lambda?
      ipt = 0; td  = 0
      do ibx1 = 1, nb
         ibx2 = dRsp_QRPA%m11%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle

         do i2=1, nd2
            ! State |2>
            ix2 = i2 + isstart(ibx2) - 1
            wf_2 = wf(:,ix2)
            xl2 = nl(ix2);  xs2 = ns(ix2)

            do i1=1, nd1
               ! State <1|
               ix1 = i1 + isstart(ibx1) - 1
               wf_1 = wf(:,ix1)
               xl1 = nl(ix1);  xs1 = ns(ix1)

               ! Matrix element
               ipt = ipt + 1

               ! Technically non-zero in td, but won't contribute to amplitude
               if (xl1 /= xl2) cycle

               ! Spin matrix element
               select case(f%K)
                  case (0)
                     if (xs1 == xs2) then
                        Sij = xs1
                     end if
                  case (1,-1)
                     if (xs1 == xs2 + 2*f%K) then
                        Sij = -f%K*sqrt(2.0_dp)
                     end if
               end select

               Rij = dRsp_QRPA%m11%elem(ipt)
               td = td + Rij*wf_1*wf_2*Sij

            end do ! i1
         end do ! i2
      end do ! ibx

      ! Get densities and NM density function
      call hfb_density_coord(IT_NEUTRON,rhon)
      call hfb_density_coord(IT_PROTON, rhop)

      rho_fac = (caux*2*cd)*(rhon + rhop)
      if (f%use_2bc(2) == 2) then
         ! Symmetric nuclear matter calculation at P=0
         rho_fac = rho_fac + tbc_nmlda_da1(rhon, rhop, 0.0_dp, .true.)
      else if (f%use_2bc(2) == 3) then
         ! Asymmetric nuclear matter calculation at P=0
         rho_fac = rho_fac + tbc_nmlda_da1(rhon, rhop, 0.0_dp, .false.)
      else if (f%use_2bc(2) == 4) then
         ! DME exchange term
         rho_fac = rho_fac + dme_exc()
      else
         rho_fac = 0
      end if

      integral_td1 = sum(td)
      integral_td2 = sum(td*rho_fac)

      ! Strip integration weights
      td = td*wdcori

      write(denfile,'(a,".den")') trim(fam_output_filename)
      call openlog(denfile, .false.)
      call writelog(" QRPA Transition Amplitude Densities")

      write(stde, '(1x,a20,es20.10e3)') "Re(EQRPA)         = ", real_eqrpa ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Im(EQRPA)         = ", imag_eqrpa ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "QRPA Norm         = ", N ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Amplitude         = ", amplitude ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Int. trans        = ", integral_td1 ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Int. trans*f(rho) = ", integral_td2 ; call writelog(stde)

      write(stde, '(6a20)') "R", "Z", "rhon", "rhop", "trans", "frho_2bc" ; call writelog(stde, .true.)
      do i1 = 1, nghl
         write(stde, '(6ES20.10E3)') 1/y(i1), z(i1), rhon(i1), rhop(i1), td(i1), rho_fac(i1)
         call writelog(stde, .true.)
      end do
      call closelog

   end subroutine transition_density

end module pnfam_solver
