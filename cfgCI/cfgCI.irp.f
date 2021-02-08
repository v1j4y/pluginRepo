      subroutine printMatrix(mat, rows, cols)
      implicit none
      BEGIN_DOC
      ! Print a 2D matrix
      END_DOC
      integer i,j
      integer,intent(in) :: rows
      integer,intent(in) :: cols
      real*8,dimension(:,:),intent(in) :: mat(rows,cols)
      print *,""
      do i=1,rows
         do j=1,cols
            write(*,'(F6.4,2X)',advance="no") mat(i,j)
         end do
         print *,""
      end do
      end subroutine

      program cfgCI
      use cfunctions
      implicit none
      BEGIN_DOC
!     TODO : Put the documentation of the program here
      END_DOC
      character*32 cc
      integer i, j, k
      integer orbp, orbq
      integer rows
      integer cols
      integer*8 MS
      print *, N_int
      !cc = "Hello"
      !print *,"Ndet=",psi_det_size, "Nconfig=",N_configuration
      !do i = 1, 6
      !   call debug_spindet(psi_configuration(N_int,1,i),N_int)
      !   call debug_spindet(psi_configuration(N_int,2,i),N_int)
      !   print *,"with mask"
      !   call debug_spindet(iand(reunion_of_act_virt_bitmask,psi_configuration(N_int,1,i)),1)
      !   call debug_spindet(iand(reunion_of_act_virt_bitmask,psi_configuration(N_int,2,i)),1)
      !end do

  integer(bit_kind) :: Icfg(N_INT,2)
  integer(bit_kind) :: alphas_Icfg(N_INT,2,200)
  integer(bit_kind) :: connectedI_alpha(N_INT,2,200)
  integer           :: idxs_connectedI_alpha(200)
  integer(bit_kind) :: psi_configuration_out(N_INT,2,400)
  real*8            :: psi_coef_out(dimBasisCSF)
  real*8            :: psi_coef_out_det(N_det)
  logical           :: psi_coef_out_init(dimBasisCSF)
  integer           :: excitationIds(2,200)
  integer           :: excitationTypes(200)
  integer  :: Nalphas_Icfg, nconnectedI, rowsikpq, colsikpq
  integer  :: extype,NSOMOalpha,Nsomoi,p,q,pmodel,qmodel
  integer :: getNSOMO
  integer :: totcolsTKI
  integer :: rowsTKI
  integer :: moi, moj, mok, mol, l,m
  real*8  :: norm_coef_cfg
  real*8  :: norm_coef_det
  real*8,dimension(:,:),allocatable :: TKI
  real*8,dimension(:,:),allocatable  :: GIJpqrs
  real*8,dimension(:,:),allocatable  :: TKIGIJ
  real*8, external :: mo_two_e_integral

  MS = 0
  norm_coef_cfg=0.d0

  psi_coef_out=0.d0
  psi_coef_out_init = .False.

  print *,"CSF basis dim=",dimBasisCSF
  do i = 1,N_configuration
     print *,i,">",psi_config_data(i,1),psi_config_data(i,2)
     call debug_spindet(psi_configuration(1,1,i),N_int)
     call debug_spindet(psi_configuration(1,2,i),N_int)
  enddo
  do i = 1,dimBasisCSF
     print *, "i=",i,"coef=",psi_coef_config(i)
     !call debug_spindet(psi_configuration(1,1,i),N_int)
     !call debug_spindet(psi_configuration(1,2,i),N_int)
     norm_coef_cfg += psi_coef_config(i)*psi_coef_config(i)
  enddo
  print *,"norm CFG = ",norm_coef_cfg
  call convertWFfromCSFtoDET(psi_coef_out,psi_coef_out_det)
  norm_coef_det=0
  do i = 1,N_det
     !print *, "i=",i,"coef=",psi_coef_out_det(i)
     norm_coef_det += psi_coef_out_det(i)*psi_coef_out_det(i)
  enddo
  print *,"norm = ",norm_coef_det, " size=",N_det

  ! Loop over all selected configurations
  do i = 1,N_configuration
     Icfg(1,1) = psi_configuration(1,1,i)
     Icfg(1,2) = psi_configuration(1,2,i)
     ! Returns all unique (checking the past) singly excited cfgs connected to I
     call obtain_associated_alphaI(i, Icfg, alphas_Icfg, Nalphas_Icfg)
     ! TODO : remove doubly excited for return
     print *,i,"Nalphas = ",Nalphas_Icfg
     call debug_spindet(Icfg(1,1),N_int)
     call debug_spindet(Icfg(1,2),N_int)
     ! Here we do 2x the loop. One to count for the size of the matrix, then we compute.
     do k = 1,Nalphas_Icfg
        print *,"Kalpha=",k
        call debug_spindet(alphas_Icfg(1,1,k),N_int)
        call debug_spindet(alphas_Icfg(1,2,k),N_int)
        ! Now generate all singly excited with respect to a given alpha CFG
        call obtain_connected_I_foralpha(i,alphas_Icfg(:,:,k),connectedI_alpha,idxs_connectedI_alpha,nconnectedI,excitationIds,excitationTypes)

        print *,k,"----> nconnected = ",nconnectedI
        totcolsTKI = 0
        rowsTKI = -1
        do j = 1,nconnectedI
           NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
           NSOMOI = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,2)
           totcolsTKI += colsikpq
           if(rowsTKI .LT. rowsikpq .AND. rowsTKI .NE. -1) then
              print *,">",j,"Something is wrong in sigma-vector", rowsTKI, rowsikpq, "(p,q)=",pmodel,qmodel,"ex=",extype,"na=",NSOMOalpha," nI=",NSOMOI
              !rowsTKI = rowsikpq
           else
              rowsTKI = rowsikpq
           endif
           !print *,"----------------alpha------"
           print *,k, Nalphas_Icfg, "idxI=",idxs_connectedI_alpha(j)
           !call debug_spindet(alphas_Icfg(1,1,k),N_int)
           !call debug_spindet(alphas_Icfg(1,2,k),N_int)
           !print *,"----------------Icfg------- Isingle=",j
           !call debug_spindet(connectedI_alpha(1,1,j),N_int)
           !call debug_spindet(connectedI_alpha(1,2,j),N_int)
           !print *,"----------------",NSOMOalpha,NSOMOI,"ex=",extype,pmodel,qmodel,"(",rowsikpq,colsikpq,")"
        end do

        !print *,"total columnTKI=",totcolsTKI
        !print *,"total rowsTKI=",rowsTKI
        ! allocate memory for table
        ! for 1 root
        ! for n roots dims = (rowsTKI,nroots,totcolsTKI)
        allocate(TKI(rowsTKI,totcolsTKI)) ! coefficients of CSF
        ! Initialize the inegral container
        ! dims : (totcolsTKI, nconnectedI)
        allocate(GIJpqrs(totcolsTKI,nconnectedI))  ! gpqrs
        allocate(TKIGIJ(rowsTKI,nconnectedI))  ! gpqrs

        TKI = 0.d0
        GIJpqrs = 0.d0
        TKIGIJ = 0.d0


        totcolsTKI = 0
        do j = 1,nconnectedI
           NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
           NSOMOI = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           !print *,j,"calling to modelspaace pq=",p,q
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
           !print *,"det a"
           !call debug_spindet(alphas_Icfg(:,1,k),1)
           !call debug_spindet(alphas_Icfg(:,2,k),1)
           !print *,"det I"
           !call debug_spindet(connectedI_alpha(:,1,j),1)
           !call debug_spindet(connectedI_alpha(:,2,j),1)
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,2)
           !print *,"j=",j,">",rowsikpq,colsikpq,"ex=",extype,"pmod(p)=",p,"qmod(q)=",q," somoI=",NSOMOI," somoa=",NSOMOalpha
           do l = 1,rowsTKI
              do m = 1,colsikpq
                 TKI(l,totcolsTKI+m) = AIJpqContainer(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,l,m) * psi_coef_config(idxs_connectedI_alpha(j)+m-1)
              enddo
           enddo
           do m = 1,colsikpq
              do l = 1,nconnectedI
                 ! <ij|kl> = (ik|jl)
                 moi = excitationIds(1,j)
                 mok = excitationIds(2,j)
                 moj = excitationIds(1,l)
                 mol = excitationIds(2,l)
                 GIJpqrs(totcolsTKI+m,l) = mo_two_e_integral(moi,moj,mok,mol)
              enddo
           enddo
           totcolsTKI += colsikpq
        end do


        print *,"TKI matrix"
        call printMatrix(TKI,rowsTKI,totcolsTKI)
        print *,"GIJpqrs matrix"
        call printMatrix(GIJpqrs,totcolsTKI,nconnectedI)

        ! Do big BLAS
        ! TODO TKI, size(TKI,1)*size(TKI,2)
        call dgemm('N','N', rowsTKI, nconnectedI, totcolsTKI, 1.d0,  &
          TKI, size(TKI,1), GIJpqrs, size(GIJpqrs,1), 0.d0, &
          TKIGIJ , size(TKIGIJ,1) )

        print *,"TKIGIJ matrix"
        call printMatrix(TKIGIJ,rowsTKI,nconnectedI)

        ! Collect the result
        totcolsTKI = 0
        do j = 1,nconnectedI
           NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
           NSOMOI     = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,2)
           !print *,">j=",j,rowsikpq,colsikpq, ">>",totcolsTKI,",",idxs_connectedI_alpha(j)
           do m = 1,colsikpq
              do l = 1,rowsTKI
                 psi_coef_out(idxs_connectedI_alpha(j)+m-1) += AIJpqContainer(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,l,m) * TKIGIJ(l,j)
                 psi_coef_out_init(idxs_connectedI_alpha(j)+m-1) = .True.
              enddo
           enddo
           totcolsTKI += colsikpq
        enddo

        deallocate(TKI) ! coefficients of CSF
        ! Initialize the inegral container
        ! dims : (totcolsTKI, nconnectedI)
        deallocate(GIJpqrs)  ! gpqrs
        deallocate(TKIGIJ)  ! gpqrs

     end do
  end do

  do i = 1,dimBasisCSF
     print *, "i=",i,"coef=",psi_coef_config(i),psi_coef_out(i)," ini?=",psi_coef_out_init(i)
  enddo

  end
