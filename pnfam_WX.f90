module pnfam_WX                      
   integer,  parameter, private :: dp = kind(1d0), pr = kind(1d0), ipr = kind(1)
contains

  subroutine effective_WX(nxy,ReX, ImX ,ReY,ImY)                 
      implicit none
      integer :: i,j  
      integer :: nxy

          
      real(dp), dimension(nxy) :: ReX, ImX, ReY, ImY     
          
      write(*,*) nxy , size(ReX)                   
  end subroutine effective_WX               
  !=======================================
  subroutine WX_cal(nxy, Ph_nb, Ph_db, Ph_readin_select, &
                    Ph_sum_db, Ph_sum_db_colm1, Ph_NumberofX, WX_sum_nd,& 
                    ReXdp, ImXdp, ReYdp, ImYdp, WX_matpnbasisN11_E101, WX_matpnbasisP11_E101, &
                    Ph_ReWX, Ph_ImWX, Ph_ReWY, Ph_ImWY, &
                    WgreenX, WgreenY, &
                    WX_iHFB_NB,WX_pnbasis_select,X_select,Y_select)       
      implicit none
      !----------------------------------------------
      ! input and output arguments
      !---------------------------------------------
      integer :: i,nxy, Ph_nb, Ph_sum_db, Ph_NumberofX, WX_sum_nd
      integer, dimension(Ph_nb, Ph_nb) :: Ph_readin_select, X_select, Y_select     
      integer, dimension(Ph_nb) :: Ph_sum_db_colm1, Ph_db         !!!Q: dimension is iHFB_NB, record sum of nd of first column minus 1 blocks.
      real(dp), dimension(nxy) :: ReXdp, ImXdp, ReYdp, ImYdp
      real, dimension(nxy) :: ReX, ImX, ReY, ImY
      Complex, dimension(2*WX_sum_nd, 2*WX_sum_nd) :: WX_matpnbasisN11_E101, WX_matpnbasisP11_E101
      real, dimension(Ph_NumberofX),intent(out) :: Ph_ReWX, Ph_ImWX, Ph_ReWY, Ph_ImWY
      !-----------------------------------------------
      ! process arguments
      !-----------------------------------------------
      !---------------------------------------------
      ! matrix PX, PPX,  
      !---------------------------------------------
      !interm variable matrix, used for both WX and WY
      !------------------------------------------
      Complex, dimension(Ph_sum_db,Ph_sum_db) :: Ph_matPdX, Ph_matPdXg, Ph_matXNc, Ph_matXNcg  
      Complex, dimension(Ph_sum_db,Ph_sum_db) :: Ph_matPdY, Ph_matPdYg, Ph_matYNc, Ph_matYNcg  
      !----------------------------------------------
      ! matrix distinguished.
      !---------------------------------------------
      Complex, dimension(Ph_sum_db,Ph_sum_db) :: Ph_matX, Ph_matY, Ph_matWX, Ph_matWY                             
      Complex, dimension(Ph_sum_db,Ph_sum_db) :: Ph_matPPdXg, Ph_matXNcgN, Ph_matPXNcg, Ph_matPdXgN   
      Complex, dimension(Ph_sum_db,Ph_sum_db) :: Ph_matPPdYg, Ph_matYNcgN, Ph_matPYNcg, Ph_matPdYgN    
      !---------------------------------------------
      complex, dimension(Ph_sum_db**2)  :: WgreenX                                                         
      complex, dimension(Ph_sum_db**2)  :: WgreenY                                                         


      !---------------------------------------------
      !interm variable matrix, real part
      !------------------------------------------
    ! real(dp), dimension(Ph_sum_db,Ph_sum_db) :: RePh_matPdX, RePh_matPdXg, RePh_matXNc, RePh_matXNcg  
    ! real(dp), dimension(Ph_sum_db,Ph_sum_db) :: RePh_matPdY, RePh_matPdYg, RePh_matYNc, RePh_matYNcg  
      !---------------------------------------------
    ! real(dp), dimension(Ph_sum_db,Ph_sum_db) :: RePh_matX, RePh_matY, RePh_matWX, RePh_matWY                             
    ! real(dp), dimension(Ph_sum_db,Ph_sum_db) :: RePh_matPPdXg, RePh_matXNcgN, RePh_matPXNcg, RePh_matPdXgN   
    ! real(dp), dimension(Ph_sum_db,Ph_sum_db) :: RePh_matPPdYg, RePh_matYNcgN, RePh_matPYNcg, RePh_matPdYgN    
    ! !---------------------------------------------
    ! real(dp), dimension(Ph_sum_db**2)  :: ReWgreenX                                                         
    ! real(dp), dimension(Ph_sum_db**2)  :: ReWgreenY                                                         
    ! !---------------------------------------------
      !interm variable matrix, imag part
      !------------------------------------------
    ! real(dp), dimension(Ph_sum_db,Ph_sum_db) :: ImPh_matPdX, ImPh_matPdXg, ImPh_matXNc, ImPh_matXNcg  
    ! real(dp), dimension(Ph_sum_db,Ph_sum_db) :: ImPh_matPdY, ImPh_matPdYg, ImPh_matYNc, ImPh_matYNcg  
    ! !---------------------------------------------
    ! real(dp), dimension(Ph_sum_db,Ph_sum_db) :: ImPh_matX, ImPh_matY, ImPh_matWX, ImPh_matWY                             
    ! real(dp), dimension(Ph_sum_db,Ph_sum_db) :: ImPh_matPPdXg, ImPh_matXNcgN, ImPh_matPXNcg, ImPh_matPdXgN   
    ! real(dp), dimension(Ph_sum_db,Ph_sum_db) :: ImPh_matPPdYg, ImPh_matYNcgN, ImPh_matPYNcg, ImPh_matPdYgN    
    ! !---------------------------------------------
    ! real(dp), dimension(Ph_sum_db**2)  :: ImWgreenX                                                         
    ! real(dp), dimension(Ph_sum_db**2)  :: ImWgreenY                                                         
      !--------------------------------------- 

      !--------------------------------------- 
      ! variable used by both matX, matY and matN11, matP11
      !---------------------------------------
      integer :: sum_block
      integer :: qqi, qqj, qqk, qql, Inselect,ct1, ct2, ct3      !!!Q: ct means count, for any
      integer :: inbcc                        !!!Q: inside block column count, used to decide which column to be assigned.
      integer :: inbrc                        !!!Q: inside block row count, used to decide whcih row element to be assigned.
      integer :: inbrr                        !!!Q: used to record how many elements already assigned to a row inside a block.
      complex,parameter :: cunit=cmplx(0.0,1.0)
      complex,parameter :: runit=cmplx(1.0,0.0)
      !---------------------------
      integer :: WX_iHFB_NB 
      integer, dimension(2*WX_iHFB_NB, 2*WX_iHFB_NB) :: WX_pnbasis_select
      integer, dimension(2*WX_iHFB_NB, 2*WX_iHFB_NB) :: AF_select, FA_select       !A ~ N or Z, F ~ X or Y
      integer, dimension(2*WX_iHFB_NB, 2*WX_iHFB_NB) :: AAF_select, AFA_select, FAA_select


      integer, dimension(2*WX_iHFB_NB, 2*WX_iHFB_NB) :: AX_select, XA_select       !A ~ N or Z, F ~ X or Y
      integer, dimension(2*WX_iHFB_NB, 2*WX_iHFB_NB) :: AAX_select, AXA_select, XAA_select
      integer, dimension(2*WX_iHFB_NB, 2*WX_iHFB_NB) :: AY_select, YA_select       !A ~ N or Z, F ~ X or Y
      integer, dimension(2*WX_iHFB_NB, 2*WX_iHFB_NB) :: AAY_select, AYA_select, YAA_select


      ReX = ReXdp 
      ImX = ImXdp
      ReY = ReYdp 
      ImY = ImYdp

      !----------------------------------------stgreen001
    ! ReWgreenX = real(WgreenX)
    ! ReWgreenY = real(WgreenY)
    ! ImWgreenX = aimag(WgreenX)
    ! ImWgreenY = aimag(WgreenY)
      !----------------------------------------endgreen001

      AF_select = matmul(WX_pnbasis_select,Ph_readin_select)
      FA_select = matmul(Ph_readin_select,WX_pnbasis_select)
      AAF_select = matmul(WX_pnbasis_select,AF_select)                   
      AFA_select = matmul(AF_select,WX_pnbasis_select)                   
      FAA_select = matmul(FA_select,WX_pnbasis_select)

      AX_select = AF_select 
      XA_select = FA_select
      AAX_select = AAF_select                    
      AXA_select = AFA_select                   
      XAA_select = FAA_select

      AY_select = AF_select 
      YA_select = FA_select
      AAY_select = AAF_select                    
      AYA_select = AFA_select                   
      YAA_select = FAA_select



     !AX_select = matmul(WX_pnbasis_select,X_select)
     !XA_select = matmul(X_select,WX_pnbasis_select)
     !AAX_select = matmul(WX_pnbasis_select,AX_select)                   
     !AXA_select = matmul(AX_select,WX_pnbasis_select)                   
     !XAA_select = matmul(XA_select,WX_pnbasis_select)


     !AY_select = AX_select 
     !YA_select = XA_select
     !AAY_select = AAX_select                    
     !AYA_select = AXA_select                   
     !YAA_select = XAA_select

 
 
 
	  !--------------------------------------
	  Ph_matPdX=0
	  Ph_matPdXg=0
	  Ph_matXNc=0
	  Ph_matXNcg=0
          !---------------
	  Ph_matPdY=0
	  Ph_matPdYg=0
	  Ph_matYNc=0
	  Ph_matYNcg=0

          !-------------------------------------
          ! for WX
          Ph_matX=0
	  Ph_matPPdXg=0
	  Ph_matXNcgN=0
	  Ph_matPXNcg=0
	  Ph_matPdXgN=0
	  Ph_matWX=0
          !-------------------------------------
          ! for WY
          Ph_matY=0
	  Ph_matPPdYg=0
	  Ph_matYNcgN=0
	  Ph_matPYNcg=0
	  Ph_matPdYgN=0
	  Ph_matWY=0
    !     !----------------------------------QQend      modify32
    !     ! real part of variables 
	! !--------------------------------------
	! RePh_matPdX=0
	! RePh_matPdXg=0
	! RePh_matXNc=0
	! RePh_matXNcg=0
    !     !---------------
	! RePh_matPdY=0
	! RePh_matPdYg=0
	! RePh_matYNc=0
	! RePh_matYNcg=0

    !     !-------------------------------------
    !     ! for WX
    !     RePh_matX=0
	! RePh_matPPdXg=0
	! RePh_matXNcgN=0
	! RePh_matPXNcg=0
	! RePh_matPdXgN=0
	! RePh_matWX=0
    !     !-------------------------------------
    !     ! for WY
    !     RePh_matY=0
	! RePh_matPPdYg=0
	! RePh_matYNcgN=0
	! RePh_matPYNcg=0
	! RePh_matPdYgN=0
	! RePh_matWY=0
    !     !----------------------------------QQend      modify32
    !   
    !     !----------------------------------QQend      modify32
    !     ! imag part of variables 
	! !--------------------------------------
	! ImPh_matPdX=0
	! ImPh_matPdXg=0
	! ImPh_matXNc=0
	! ImPh_matXNcg=0
    !     !---------------
	! ImPh_matPdY=0
	! ImPh_matPdYg=0
	! ImPh_matYNc=0
	! ImPh_matYNcg=0

    !     !-------------------------------------
    !     ! for WX
    !     ImPh_matX=0
	! ImPh_matPPdXg=0
	! ImPh_matXNcgN=0
	! ImPh_matPXNcg=0
	! ImPh_matPdXgN=0
	! ImPh_matWX=0
    !     !-------------------------------------
    !     ! for WY
    !     ImPh_matY=0
	! ImPh_matPPdYg=0
	! ImPh_matYNcgN=0
	! ImPh_matPYNcg=0
	! ImPh_matPdYgN=0
	! ImPh_matWY=0
    !     !----------------------------------QQend      modify32
!

          !-------------------------
          ! asssign array X to matX
          !-------------------------
          !
          ct1=1
          ct3=0
          inbcc=0
          inbrc=0
          inbrr=0
          do qqj=1,Ph_nb
             do qqi=1,Ph_nb
                if(Ph_readin_select(qqi,qqj)==1)  then
               !if(X_select(qqi,qqj)==1)  then
                   inbcc=Ph_sum_db_colm1(qqj)+1
                   inbrc=Ph_sum_db_colm1(qqi)
                   ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
                      do ct2=ct1, ct3 
                          inbrc=inbrc+1
                          inbrr=inbrr+1
                          Ph_matX(inbrc, inbcc)=ReX(ct2)+cunit*ImX(ct2)             
                          Ph_matY(inbrc, inbcc)=ReY(ct2)+cunit*ImY(ct2)       !modify35          
                        if(inbrr.eq. Ph_db(qqi)) then
                           inbcc=inbcc+1        !!!Q: col increase by 1
                           inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
                           inbrr=0              !!!Q: reset inbrr
                        end if
                      end do
                      ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
                end if
                !rite(*,*) 'assign to matX, matY', qqj, qqi
             end do
          end do


        ! ct1=1
        ! ct3=0
        ! inbcc=0
        ! inbrc=0
        ! inbrr=0
        ! do qqj=1,Ph_nb
        !    do qqi=1,Ph_nb
        !       if(Ph_readin_select(qqi,qqj)==1)  then
        !      !if(Y_select(qqi,qqj)==1)  then
        !          inbcc=Ph_sum_db_colm1(qqj)+1
        !          inbrc=Ph_sum_db_colm1(qqi)
        !          ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
        !             do ct2=ct1, ct3 
        !                 inbrc=inbrc+1
        !                 inbrr=inbrr+1
        !                !Ph_matX(inbrc, inbcc)=ReX(ct2)+cunit*ImX(ct2)             
        !                 Ph_matY(inbrc, inbcc)=ReY(ct2)+cunit*ImY(ct2)       !modify35          
        !               if(inbrr.eq. Ph_db(qqi)) then
        !                  inbcc=inbcc+1        !!!Q: col increase by 1
        !                  inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
        !                  inbrr=0              !!!Q: reset inbrr
        !               end if
        !             end do
        !             ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
        !       end if
        !       !rite(*,*) 'assign to matX, matY', qqj, qqi
        !    end do
        ! end do


          !-----------------------------------------------------------------------------E101     modify35
          ! X,Y part
          !-----------------------------------------------QQst          modify25
          ! X4.1.1 construct Ph_matPdF
          !-------------------------------
          Ph_matPdX=0 
          Ph_matPdY=0 
          Ph_matPdXg=0 
          Ph_matPdYg=0 
          !
          ct1=1
          ct3=0
          inbcc=0
          inbrc=0
          inbrr=0
          do qqj=1,Ph_nb
             do qqi=1,Ph_nb
                if(AX_select(qqi,qqj).gt.0)  then

                   ! choose the nozero block of X.
                   do qql = 1, Ph_nb
                      if (Ph_readin_select(qql, qqj) .Gt. 0) then 

                   inbcc=Ph_sum_db_colm1(qqj)+1
                   inbrc=Ph_sum_db_colm1(qqi)
                   ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
                      do ct2=ct1, ct3 
                          inbrc=inbrc+1
                          inbrr=inbrr+1
                        ! do qqk=1,Ph_sum_db
                          !This should be the nozero blocks of X.
                          do qqk=Ph_sum_db_colm1(qql)+1,& 
                                 Ph_sum_db_colm1(qql)+ Ph_db(qql)
                             Ph_matPdX(inbrc,inbcc)= Ph_matPdX(inbrc,inbcc) + &
                                                       conjg(WX_matpnbasisP11_E101(qqk,inbrc))*Ph_matX(qqk,inbcc) 
                             Ph_matPdY(inbrc,inbcc)= Ph_matPdY(inbrc,inbcc) + &
                                                       WX_matpnbasisP11_E101(qqk,inbrc)*Ph_matY(qqk,inbcc)

                          end do
                        if(inbrr.eq. Ph_db(qqi)) then
                           inbcc=inbcc+1        !!!Q: col increase by 1
                           inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
                           inbrr=0              !!!Q: reset inbrr
                        end if
                      end do
                      ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
                      end if !(Ph_readin_select(qql, qqj) .eq.1)
                   end do !qql
                end if
                !rite(*,*)'4.1.1', qqj, qqi
             end do
          end do
          !
        ! ct1=1
        ! ct3=0
        ! inbcc=0
        ! inbrc=0
        ! inbrr=0
        ! do qqj=1,Ph_nb
        !    do qqi=1,Ph_nb
        !       if(AY_select(qqi,qqj).gt.0)  then
        !          inbcc=Ph_sum_db_colm1(qqj)+1
        !          inbrc=Ph_sum_db_colm1(qqi)
        !          ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
        !             do ct2=ct1, ct3 
        !                 inbrc=inbrc+1
        !                 inbrr=inbrr+1
        !               ! Ph_matX(inbrc, inbcc)=ReX(ct2)+cunit*ImX(ct2)             
        !               ! Ph_matY(inbrc, inbcc)=ReY(ct2)+cunit*ImY(ct2)       !modify35          
        !                 do qqk=1,Ph_sum_db
        !                   !Ph_matPdX(inbrc,inbcc)= Ph_matPdX(inbrc,inbcc) + &
        !                   !                          WX_matpnbasisP11_E101(qqk,inbrc)*Ph_matX(qqk,inbcc) 
        !                    Ph_matPdY(inbrc,inbcc)= Ph_matPdY(inbrc,inbcc) + &
        !                                              conjg(WX_matpnbasisP11_E101(qqk,inbrc))*conjg(Ph_matY(qqk,inbcc))

        !                 end do
        !               if(inbrr.eq. Ph_db(qqi)) then
        !                  inbcc=inbcc+1        !!!Q: col increase by 1
        !                  inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
        !                  inbrr=0              !!!Q: reset inbrr
        !               end if
        !             end do
        !             ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
        !       end if
        !       !rite(*,*)'4.1.1', qqj, qqi
        !    end do
        ! end do



          !-------------------------------
          ! X4.1.2 construct Ph_matPdFg from Ph_matPdF
          !-------------------------------
          i = 0
          do qqj=1,Ph_sum_db
             do qqi=1,Ph_sum_db
                i = i +1
                Ph_matPdXg(qqi,qqj)=Ph_matPdX(qqi,qqj)*WgreenX(i)
                Ph_matPdYg(qqi,qqj)=Ph_matPdY(qqi,qqj)*WgreenY(i)
             end do
             !rite(*,*) '4.1.2', qqj
          end do
          !-------------------------------
          ! X4.1.3 constrcut Ph_matPPdXg from Ph_matPdFg
          !-------------------------------

        ! ct1=1
        ! ct3=0
        ! inbcc=0
        ! inbrc=0
        ! inbrr=0
        ! do qqj=1,Ph_nb
        !    do qqi=1,Ph_nb
        !      !if(Ph_readin_select(qqi,qqj).gt.0)  then
        !       if(AAX_select(qqi,qqj).gt.0)  then
        !          inbcc=Ph_sum_db_colm1(qqj)+1
        !          inbrc=Ph_sum_db_colm1(qqi)
        !          ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
        !             do ct2=ct1, ct3 
        !                 inbrc=inbrc+1
        !                 inbrr=inbrr+1
        !                 do qqk=1,Ph_sum_db
        !                    Ph_matPPdXg(inbrc,inbcc)= Ph_matPPdXg(inbrc,inbcc) + &
        !                                              WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matPdXg(qqk,inbcc) 
        !                    Ph_matPPdYg(inbrc,inbcc)= Ph_matPPdYg(inbrc,inbcc) + &
        !                                              WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matPdYg(qqk,inbcc) 
        !                 end do
        !               if(inbrr.eq. Ph_db(qqi)) then
        !                  inbcc=inbcc+1        !!!Q: col increase by 1
        !                  inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
        !                  inbrr=0              !!!Q: reset inbrr
        !               end if
        !             end do
        !             ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
        !       end if
        !       !rite(*,*) '4.1.3', qqj, qqi           
        !    end do
        ! end do
          !
        ! ct1=1
        ! ct3=0
        ! inbcc=0
        ! inbrc=0
        ! inbrr=0
        ! do qqj=1,Ph_nb
        !    do qqi=1,Ph_nb
        !      !if(Ph_readin_select(qqi,qqj).gt.0)  then
        !       if(AAY_select(qqi,qqj).gt.0)  then
        !          inbcc=Ph_sum_db_colm1(qqj)+1
        !          inbrc=Ph_sum_db_colm1(qqi)
        !          ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
        !             do ct2=ct1, ct3 
        !                 inbrc=inbrc+1
        !                 inbrr=inbrr+1
        !                 do qqk=1,Ph_sum_db
        !                   !Ph_matPPdXg(inbrc,inbcc)= Ph_matPPdXg(inbrc,inbcc) + &
        !                   !                          WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matPdXg(qqk,inbcc) 
        !                    Ph_matPPdYg(inbrc,inbcc)= Ph_matPPdYg(inbrc,inbcc) + &
        !                                              WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matPdYg(qqk,inbcc) 
        !                 end do
        !               if(inbrr.eq. Ph_db(qqi)) then
        !                  inbcc=inbcc+1        !!!Q: col increase by 1
        !                  inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
        !                  inbrr=0              !!!Q: reset inbrr
        !               end if
        !             end do
        !             ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
        !       end if
        !       !rite(*,*) '4.1.3', qqj, qqi           
        !    end do
        ! end do



          !
          !-------------------------------
          ! X4.2.1 construct Ph_matFNc
          !-------------------------------
          !
          ct1=1
          ct3=0
          inbcc=0
          inbrc=0
          inbrr=0
          do qqj=1,Ph_nb
             do qqi=1,Ph_nb
                if(XA_select(qqi,qqj).gt.0)  then

                   ! choose the nozero block of X.
                   do qql = 1, Ph_nb
                      if (Ph_readin_select(qqi, qql) .Gt. 0) then 


                   inbcc=Ph_sum_db_colm1(qqj)+1
                   inbrc=Ph_sum_db_colm1(qqi)
                   ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
                      do ct2=ct1, ct3 
                          inbrc=inbrc+1
                          inbrr=inbrr+1
                         !do qqk=1,Ph_sum_db
                          do qqk=Ph_sum_db_colm1(qql)+1, & 
                                 Ph_sum_db_colm1(qql)+Ph_db(qql)
                             Ph_matXNc(inbrc,inbcc)= Ph_matXNc(inbrc,inbcc) + &
                                                     Ph_matX(inbrc,qqk)*conjg(WX_matpnbasisN11_E101(qqk,inbcc)) 
                             Ph_matYNc(inbrc,inbcc)= Ph_matYNc(inbrc,inbcc) + &
                                                     Ph_matY(inbrc,qqk)*WX_matpnbasisN11_E101(qqk,inbcc) 

                          end do
                        if(inbrr.eq. Ph_db(qqi)) then
                           inbcc=inbcc+1        !!!Q: col increase by 1
                           inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
                           inbrr=0              !!!Q: reset inbrr
                        end if
                      end do
                      ct1=ct1+Ph_db(qqj)*Ph_db(qqi)

                      end if !(Ph_readin_select(qqi, qql) .eq.1)
                   end do !qql
 
                end if
                !rite(*,*) '4.2.1',qqj, qqi
             end do
          end do
          !
       !  ct1=1
       !  ct3=0
       !  inbcc=0
       !  inbrc=0
       !  inbrr=0
       !  do qqj=1,Ph_nb
       !     do qqi=1,Ph_nb
       !        if(YA_select(qqi,qqj).gt.0)  then
       !           inbcc=Ph_sum_db_colm1(qqj)+1
       !           inbrc=Ph_sum_db_colm1(qqi)
       !           ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
       !              do ct2=ct1, ct3 
       !                  inbrc=inbrc+1
       !                  inbrr=inbrr+1
       !                  do qqk=1,Ph_sum_db
       !                    !Ph_matXNc(inbrc,inbcc)= Ph_matXNc(inbrc,inbcc) + &
       !                    !                        Ph_matX(inbrc,qqk)*WX_matpnbasisN11_E101(qqk,inbcc) 
       !                     Ph_matYNc(inbrc,inbcc)= Ph_matYNc(inbrc,inbcc) + &
       !                                             conjg(Ph_matY(inbrc,qqk))*conjg(WX_matpnbasisN11_E101(qqk,inbcc)) 

       !                  end do
       !                if(inbrr.eq. Ph_db(qqi)) then
       !                   inbcc=inbcc+1        !!!Q: col increase by 1
       !                   inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
       !                   inbrr=0              !!!Q: reset inbrr
       !                end if
       !              end do
       !              ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
       !        end if
       !        !rite(*,*) '4.2.1',qqj, qqi
       !     end do
       !  end do



          !-------------------------------
          ! X4.2.2 construct Ph_matFNcg
          !-------------------------------
          i = 0
          do qqj=1,Ph_sum_db
             do qqi=1,Ph_sum_db
                i = i+1
                Ph_matXNcg(qqi,qqj)=Ph_matXNc(qqi,qqj)*WgreenX(i)
                Ph_matYNcg(qqi,qqj)=Ph_matYNc(qqi,qqj)*WgreenY(i)
             end do
             !rite(*,*) '4.2.2',qqj
          end do
          !      
          !-------------------------------
          ! X4.2.3 construct Ph_matXNcgN
          !-------------------------------
        ! ct1=1
        ! ct3=0
        ! inbcc=0
        ! inbrc=0
        ! inbrr=0
        ! do qqj=1,Ph_nb
        !    do qqi=1,Ph_nb
        !      !if(Ph_readin_select(qqi,qqj)==1)  then
        !       if(XAA_select(qqi,qqj).gt.0)  then
        !          inbcc=Ph_sum_db_colm1(qqj)+1
        !          inbrc=Ph_sum_db_colm1(qqi)
        !          ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
        !             do ct2=ct1, ct3 
        !                 inbrc=inbrc+1
        !                 inbrr=inbrr+1
        !                 do qqk=1,Ph_sum_db
        !                    Ph_matXNcgN(inbrc,inbcc)= Ph_matXNcgN(inbrc,inbcc) + &
        !                                              Ph_matXNcg(inbrc,qqk)*WX_matpnbasisN11_E101(inbcc,qqk)    
        !                    Ph_matYNcgN(inbrc,inbcc)= Ph_matYNcgN(inbrc,inbcc) + &
        !                                              Ph_matYNcg(inbrc,qqk)*WX_matpnbasisN11_E101(inbcc,qqk)  
        !                 end do
        !               if(inbrr.eq. Ph_db(qqi)) then
        !                  inbcc=inbcc+1        !!!Q: col increase by 1
        !                  inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
        !                  inbrr=0              !!!Q: reset inbrr
        !               end if
        !             end do
        !             ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
        !       end if
        !       !rite(*,*) '4.2.3', qqj, qqi
        !    end do
        ! end do
          !
       !  ct1=1
       !  ct3=0
       !  inbcc=0
       !  inbrc=0
       !  inbrr=0
       !  do qqj=1,Ph_nb
       !     do qqi=1,Ph_nb
       !       !if(Ph_readin_select(qqi,qqj)==1)  then
       !        if(YAA_select(qqi,qqj).gt.0)  then
       !           inbcc=Ph_sum_db_colm1(qqj)+1
       !           inbrc=Ph_sum_db_colm1(qqi)
       !           ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
       !              do ct2=ct1, ct3 
       !                  inbrc=inbrc+1
       !                  inbrr=inbrr+1
       !                  do qqk=1,Ph_sum_db
       !                    !Ph_matXNcgN(inbrc,inbcc)= Ph_matXNcgN(inbrc,inbcc) + &
       !                    !                          Ph_matXNcg(inbrc,qqk)*WX_matpnbasisN11_E101(inbcc,qqk)    
       !                     Ph_matYNcgN(inbrc,inbcc)= Ph_matYNcgN(inbrc,inbcc) + &
       !                                               Ph_matYNcg(inbrc,qqk)*WX_matpnbasisN11_E101(inbcc,qqk)  
       !                  end do
       !                if(inbrr.eq. Ph_db(qqi)) then
       !                   inbcc=inbcc+1        !!!Q: col increase by 1
       !                   inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
       !                   inbrr=0              !!!Q: reset inbrr
       !                end if
       !              end do
       !              ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
       !        end if
       !        !rite(*,*) '4.2.3', qqj, qqi
       !     end do
       !  end do



          !
          !-------------------------------
          ! X4.3.1 construct Ph_matFNc 
          !-------------------------------
          ! Ph_matXNc constructed in 4.2.1
          !-------------------------------
          ! X4.3.2 construct Ph_matFNcg
          !-------------------------------
          ! Ph_matXNcg constructed in 4.2.2
          !-------------------------------
          ! X4.3.3 construct Ph_matPXNcg
          !-------------------------------
        ! ct1=1
        ! ct3=0
        ! inbcc=0
        ! inbrc=0
        ! inbrr=0
        ! do qqj=1,Ph_nb
        !    do qqi=1,Ph_nb
        !      !if(Ph_readin_select(qqi,qqj)==1)  then
        !       if(AXA_select(qqi,qqj).gt.0)  then
        !          inbcc=Ph_sum_db_colm1(qqj)+1
        !          inbrc=Ph_sum_db_colm1(qqi)
        !          ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
        !             do ct2=ct1, ct3 
        !                 inbrc=inbrc+1
        !                 inbrr=inbrr+1
        !                 do qqk=1,Ph_sum_db
        !                    Ph_matPXNcg(inbrc,inbcc)= Ph_matPXNcg(inbrc,inbcc) + &
        !                                              WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matXNcg(qqk,inbcc) 
        !                    Ph_matPYNcg(inbrc,inbcc)= Ph_matPYNcg(inbrc,inbcc) + &
        !                                              WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matYNcg(qqk,inbcc) 
        !                 end do
        !               if(inbrr.eq. Ph_db(qqi)) then
        !                  inbcc=inbcc+1        !!!Q: col increase by 1
        !                  inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
        !                  inbrr=0              !!!Q: reset inbrr
        !               end if
        !             end do
        !             ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
        !       end if
        !       !rite(*,*) '4.3.3', qqj, qqi
        !    end do
        ! end do
        ! !
        ! ct1=1
        ! ct3=0
        ! inbcc=0
        ! inbrc=0
        ! inbrr=0
        ! do qqj=1,Ph_nb
        !    do qqi=1,Ph_nb
        !      !if(Ph_readin_select(qqi,qqj)==1)  then
        !       if(AYA_select(qqi,qqj).gt.0)  then
        !          inbcc=Ph_sum_db_colm1(qqj)+1
        !          inbrc=Ph_sum_db_colm1(qqi)
        !          ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
        !             do ct2=ct1, ct3 
        !                 inbrc=inbrc+1
        !                 inbrr=inbrr+1
        !                 do qqk=1,Ph_sum_db
        !                   !Ph_matPXNcg(inbrc,inbcc)= Ph_matPXNcg(inbrc,inbcc) + &
        !                   !                          WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matXNcg(qqk,inbcc) 
        !                    Ph_matPYNcg(inbrc,inbcc)= Ph_matPYNcg(inbrc,inbcc) + &
        !                                              WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matYNcg(qqk,inbcc) 
        !                 end do
        !               if(inbrr.eq. Ph_db(qqi)) then
        !                  inbcc=inbcc+1        !!!Q: col increase by 1
        !                  inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
        !                  inbrr=0              !!!Q: reset inbrr
        !               end if
        !             end do
        !             ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
        !       end if
        !       !rite(*,*) '4.3.3', qqj, qqi
        !    end do
        ! end do


          !
          !-------------------------------
          ! X4.4.1 construct Ph_matPdX 
          !-------------------------------
          ! Ph_matPdF constructed in 4.1.1
          !-------------------------------
          ! x4.4.2 construct Ph_matPdFg
          !-------------------------------
          ! Ph_matPdFg constructed in 4.1.2
          !-------------------------------
          ! X4.4.3 construct Ph_matPdXgN
          !-------------------------------
          ct1=1
          ct3=0
          inbcc=0
          inbrc=0
          inbrr=0
          do qqj=1,Ph_nb
             do qqi=1,Ph_nb
                if(Ph_readin_select(qqi,qqj)==1)  then

                   ! choose the nozero block of X.
                   do qql = 1, Ph_nb
                      if (WX_pnbasis_select(qqi, qql) .Gt. 0) then 


               !if(AXA_select(qqi,qqj).gt.0)  then
                   inbcc=Ph_sum_db_colm1(qqj)+1
                   inbrc=Ph_sum_db_colm1(qqi)
                   ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
                      do ct2=ct1, ct3 
                          inbrc=inbrc+1
                          inbrr=inbrr+1
                         !do qqk=1,Ph_sum_db
                          do qqk=Ph_sum_db_colm1(qql)+1, & 
                                 Ph_sum_db_colm1(qql)+Ph_db(qql)
 
                             Ph_matPPdXg(inbrc,inbcc)= Ph_matPPdXg(inbrc,inbcc) + &
                                                       WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matPdXg(qqk,inbcc) 
                             Ph_matPPdYg(inbrc,inbcc)= Ph_matPPdYg(inbrc,inbcc) + &
                                                       WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matPdYg(qqk,inbcc) 
 

                             Ph_matPXNcg(inbrc,inbcc)= Ph_matPXNcg(inbrc,inbcc) + &
                                                       WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matXNcg(qqk,inbcc) 
                             Ph_matPYNcg(inbrc,inbcc)= Ph_matPYNcg(inbrc,inbcc) + &
                                                       WX_matpnbasisP11_E101(inbrc,qqk)*Ph_matYNcg(qqk,inbcc) 
 
                          end do
                        if(inbrr.eq. Ph_db(qqi)) then
                           inbcc=inbcc+1        !!!Q: col increase by 1
                           inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
                           inbrr=0              !!!Q: reset inbrr
                        end if
                      end do
                      ct1=ct1+Ph_db(qqj)*Ph_db(qqi)

                      end if !(Ph_readin_select(qqi, qql) .eq.1)
                   end do !qql
 
                end if
                !rite(*,*) '4.4.3', qqj, qqi
             end do
          end do


          ct1=1
          ct3=0
          inbcc=0
          inbrc=0
          inbrr=0
          do qqj=1,Ph_nb
             do qqi=1,Ph_nb
                if(Ph_readin_select(qqi,qqj)==1)  then

                   ! choose the nozero block of X.
                   do qql = 1, Ph_nb
                      if (WX_pnbasis_select(qql, qqj) .Gt. 0) then 


               !if(AXA_select(qqi,qqj).gt.0)  then
                   inbcc=Ph_sum_db_colm1(qqj)+1
                   inbrc=Ph_sum_db_colm1(qqi)
                   ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
                      do ct2=ct1, ct3 
                          inbrc=inbrc+1
                          inbrr=inbrr+1
                         !do qqk=1,Ph_sum_db
                          do qqk=Ph_sum_db_colm1(qql)+1, & 
                                 Ph_sum_db_colm1(qql)+Ph_db(qql)
 

                             Ph_matXNcgN(inbrc,inbcc)= Ph_matXNcgN(inbrc,inbcc) + &
                                                       Ph_matXNcg(inbrc,qqk)*WX_matpnbasisN11_E101(inbcc,qqk)    
                             Ph_matYNcgN(inbrc,inbcc)= Ph_matYNcgN(inbrc,inbcc) + &
                                                       Ph_matYNcg(inbrc,qqk)*WX_matpnbasisN11_E101(inbcc,qqk)  
 

                             Ph_matPdXgN(inbrc,inbcc)= Ph_matPdXgN(inbrc,inbcc) + &
                                                       Ph_matPdXg(inbrc,qqk)*WX_matpnbasisN11_E101(inbcc,qqk)
                             Ph_matPdYgN(inbrc,inbcc)= Ph_matPdYgN(inbrc,inbcc) + &
                                                       Ph_matPdYg(inbrc,qqk)*WX_matpnbasisN11_E101(inbcc,qqk)
                          end do
                        if(inbrr.eq. Ph_db(qqi)) then
                           inbcc=inbcc+1        !!!Q: col increase by 1
                           inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
                           inbrr=0              !!!Q: reset inbrr
                        end if
                      end do
                      ct1=ct1+Ph_db(qqj)*Ph_db(qqi)

                      end if !(Ph_readin_select(qqi, qql) .eq.1)
                   end do !qql
 
                end if
                !rite(*,*) '4.4.3', qqj, qqi
             end do
          end do
 
         !
        ! ct1=1
        ! ct3=0
        ! inbcc=0
        ! inbrc=0
        ! inbrr=0
        ! do qqj=1,Ph_nb
        !    do qqi=1,Ph_nb
        !      !if(Ph_readin_select(qqi,qqj)==1)  then
        !       if(AYA_select(qqi,qqj).gt.0)  then
        !          inbcc=Ph_sum_db_colm1(qqj)+1
        !          inbrc=Ph_sum_db_colm1(qqi)
        !          ct3=ct3+Ph_db(qqj)*Ph_db(qqi)
        !             do ct2=ct1, ct3 
        !                 inbrc=inbrc+1
        !                 inbrr=inbrr+1
        !                 do qqk=1,Ph_sum_db
        !                   !Ph_matPdXgN(inbrc,inbcc)= Ph_matPdXgN(inbrc,inbcc) + &
        !                   !                          Ph_matPdXg(inbrc,qqk)*WX_matpnbasisN11_E101(inbcc,qqk)
        !                    Ph_matPdYgN(inbrc,inbcc)= Ph_matPdYgN(inbrc,inbcc) + &
        !                                              Ph_matPdYg(inbrc,qqk)*WX_matpnbasisN11_E101(inbcc,qqk)
        !                 end do
        !               if(inbrr.eq. Ph_db(qqi)) then
        !                  inbcc=inbcc+1        !!!Q: col increase by 1
        !                  inbrc=Ph_sum_db_colm1(qqi)   !!!Q: reset inbrc   
        !                  inbrr=0              !!!Q: reset inbrr
        !               end if
        !             end do
        !             ct1=ct1+Ph_db(qqj)*Ph_db(qqi)
        !       end if
        !       !rite(*,*) '4.4.3', qqj, qqi
        !    end do
        ! end do



          !-------------------------------
          Ph_matWX= Ph_matWX +Ph_matPPdXg +Ph_matXNcgN +Ph_matPXNcg + Ph_matPdXgN
          Ph_matWY= Ph_matWY +Ph_matPPdYg +Ph_matYNcgN +Ph_matPYNcg + Ph_matPdYgN
          !----------------------------------------------------------------------E101    

          !-----------------------------------------------------------------------------002     modify35

          
          !----------------------------------------------------------------------002
          !-------------------------------
 
          !-------------------------------
          ! restore Ph_matWX and ph_matWY to array form, block info is pnfam block structure.
          ! Ph_readin_select
          !-------------------------------
          inbcc=0
          inbrc=0
          ct3=0
          do qqj=1,Ph_nb      
             do qqi=1,Ph_nb
                if(Ph_readin_select(qqi,qqj)==1)  then
               !if(X_select(qqi,qqj)==1)  then
                   inbcc=Ph_sum_db_colm1(qqj)
                   inbrc=Ph_sum_db_colm1(qqi)
                   do ct1=1, Ph_db(qqj)     
                      do ct2=1, Ph_db(qqi)
                         ct3=ct3+1
                         !--------------------
                         Ph_ReWX(ct3)=real(Ph_matWX(inbrc+ct2,inbcc+ct1))
                         Ph_ImWX(ct3)=aimag(Ph_matWX(inbrc+ct2,inbcc+ct1))
                         !------------------- 
                         Ph_ReWY(ct3)=real(Ph_matWY(inbrc+ct2,inbcc+ct1))
                         Ph_ImWY(ct3)=aimag(Ph_matWY(inbrc+ct2,inbcc+ct1))

                      end do
                   end do
                end if
             end do
          end do
         !
       !  inbcc=0
       !  inbrc=0
       !  ct3=0
       !  do qqj=1,Ph_nb      
       !     do qqi=1,Ph_nb
       !        if(Ph_readin_select(qqi,qqj)==1)  then
       !       !if(Y_select(qqi,qqj)==1)  then
       !           inbcc=Ph_sum_db_colm1(qqj)
       !           inbrc=Ph_sum_db_colm1(qqi)
       !           do ct1=1, Ph_db(qqj)     
       !              do ct2=1, Ph_db(qqi)
       !                 ct3=ct3+1
       !                 !--------------------
       !                !Ph_ReWX(ct3)=real(Ph_matWX(inbrc+ct2,inbcc+ct1))
       !                !Ph_ImWX(ct3)=aimag(Ph_matWX(inbrc+ct2,inbcc+ct1))
       !                 !------------------- 
       !                 Ph_ReWY(ct3)=real(Ph_matWY(inbrc+ct2,inbcc+ct1))
       !                 Ph_ImWY(ct3)=aimag(Ph_matWY(inbrc+ct2,inbcc+ct1))

       !              end do
       !           end do
       !        end if
       !     end do
       !  end do


  end subroutine WX_cal            

  !=======================================
  subroutine WgreenX_func(Ph_sum_db, Ep, En, omega, WX_realomega, width, WgreenX) 
     implicit none
     integer :: qqi, qqj
     integer :: i,Ph_sum_db
     real(dp), dimension(Ph_sum_db) :: Ep, En
     real(dp) :: width, omega, WX_realomega
     complex, dimension(Ph_sum_db**2), intent(out)  :: WgreenX(:)                                                         
     complex,parameter :: cunit=cmplx(0.0_dp,1.0_dp,kind=dp)
     real :: denominate = 0

   ! if(WX_realomega .lt. 1.0) then 
   !    denominate = 0.1/WX_realomega
   ! end if
 
     i = 0
     do qqj=1,Ph_sum_db
         do qqi=1,Ph_sum_db
            i = i+1
          ! denominate = omega-WX_realomega-Ep(qqi)-En(qqj)
          ! if (abs(denominate) .gt.10) then
          !     WgreenX(i) = 0
          ! else 
                WgreenX(i) = 1/(omega+denominate*cunit+width*cunit-WX_realomega-Ep(qqi)-En(qqj))-1/(denominate*cunit+width*cunit-WX_realomega-Ep(qqi)-En(qqj))                                           
           !WgreenX(i) = 1/(omega+0.2*width*cunit-WX_realomega-Ep(qqi)-En(qqj))                                          
            !rite(*,*) 'WgreenX',qqj, qqi
          ! end if
         end do
     end do
     !
  end subroutine WgreenX_func
  !=====================================
  subroutine WgreenY_func(Ph_sum_db, Ep, En, omega, WX_realomega, width, WgreenY) 
     implicit none
     integer :: qqi, qqj
     integer :: i,Ph_sum_db
     real(dp), dimension(Ph_sum_db) :: Ep, En
     real(dp) :: width, omega, WX_realomega
     complex, dimension(Ph_sum_db**2), intent(out)  :: WgreenY(:)                                                         
     complex,parameter :: cunit=cmplx(0.0_dp,1.0_dp,kind=dp)
     real :: denominate = 0
 
    !if(WX_realomega .lt. 1.0) then 
    !   denominate = 0.1/WX_realomega
    !end if
 
     i = 0
     do qqj=1,Ph_sum_db
         do qqi=1,Ph_sum_db
            i = i+1
          ! denominate = omega+WX_realomega+Ep(qqi)+En(qqj)
          ! if (abs(denominate) .gt.10) then 
          !     WgreenY(i) = 0
          ! else
                WgreenY(i) = 1/(-omega-denominate*cunit-width*cunit-WX_realomega-Ep(qqi)-En(qqj))-1/(-denominate*cunit-width*cunit-WX_realomega-Ep(qqi)-En(qqj))                                          
           !WgreenY(i) = 1/(-omega-0.2*width*cunit-WX_realomega-Ep(qqi)-En(qqj))                                          
            !rite(*,*) 'WgreenY',qqj, qqi
          ! end if
         end do
     end do
     !
  end subroutine WgreenY_func


end module pnfam_WX                   
