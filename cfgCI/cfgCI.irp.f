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
      print *,N_int
      cc = "Hello"
      print *,"Ndet=",psi_det_size, "Nconfig=",N_configuration
      do i = 1, 6
         call debug_spindet(psi_configuration(N_int,1,i),N_int)
         call debug_spindet(psi_configuration(N_int,2,i),N_int)
         print *,"with mask"
         call debug_spindet(iand(reunion_of_act_virt_bitmask,psi_configuration(N_int,1,i)),1)
         call debug_spindet(iand(reunion_of_act_virt_bitmask,psi_configuration(N_int,2,i)),1)
      end do

  integer(bit_kind) :: Icfg(N_INT,2)
  integer(bit_kind) :: alphas_Icfg(N_INT,2,200)
  integer(bit_kind) :: connectedI_alpha(N_INT,2,200)
  integer(bit_kind) :: psi_configuration_out(N_INT,2,400)
  real*8            :: psi_coef_out(400)
  integer           :: excitationIds(2,200)
  integer           :: excitationTypes(200)
  integer  :: Nalphas_Icfg, nconnectedI, rowsikpq, colsikpq
  integer  :: extype,NSOMOalpha,Nsomoi,p,q,pmodel,qmodel
  integer :: getNSOMO
  integer :: totcolsTKI
  integer :: rowsTKI
  integer :: moi, moj, mok, mol, l,m
  real*8,dimension(:,:),allocatable :: TKI
  real*8,dimension(:,:),allocatable  :: GIJpqrs
  real*8,dimension(:,:),allocatable  :: TKIGIJ
  real*8, external :: mo_two_e_integral

  MS = 0

  ! Loop over all selected configurations
  do i = 1,N_configuration
     Icfg(1,1) = psi_configuration(1,1,i)
     Icfg(1,2) = psi_configuration(1,2,i)
     ! Returns all unique (checking the past) singly excited cfgs connected to I
     call obtain_associated_alphaI(i, Icfg, alphas_Icfg, Nalphas_Icfg)
     ! TODO : remove doubly excited for return
     print *,i,"Nalphas = ",Nalphas_Icfg
     ! Here we do 2x the loop. One to count for the size of the matrix, then we compute.
     do k = 1,Nalphas_Icfg
        ! Now generate all singly excited with respect to a given alpha CFG
        call obtain_connected_I_foralpha(i,alphas_Icfg(:,:,k),connectedI_alpha,nconnectedI,excitationIds,excitationTypes)

        print *,k,"----> nconnected = ",nconnectedI
        totcolsTKI = 0
        do j = 1,nconnectedI
           NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
           NSOMOI = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha+1,NSOMOI+1,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha+1,NSOMOI+1,extype,pmodel,qmodel,2)
           totcolsTKI += colsikpq
           if(rowsTKI .LT. rowsikpq) rowsTKI = rowsikpq
           !print *,"----------------alpha------"
           !print *,k, Nalphas_Icfg
           !call debug_spindet(alphas_Icfg(1,1,k),N_int)
           !call debug_spindet(alphas_Icfg(1,2,k),N_int)
           !print *,"----------------Icfg------- Isingle=",j
           !call debug_spindet(connectedI_alpha(1,1,j),N_int)
           !call debug_spindet(connectedI_alpha(1,2,j),N_int)
           print *,"----------------",NSOMOalpha,NSOMOI,pmodel,qmodel,"(",rowsikpq,colsikpq,")"
        end do

        print *,"total columns=",totcolsTKI
        ! allocate memory for table
        ! for 1 root
        ! for n roots dims = (rowsTKI,nroots,totcolsTKI)
        allocate(TKI(rowsTKI,totcolsTKI)) ! coefficients of CSF
        ! Initialize the inegral container
        ! dims : (totcolsTKI, nconnectedI)
        allocate(GIJpqrs(totcolsTKI,nconnectedI))  ! gpqrs
        allocate(TKIGIJ(rowsTKI,nconnectedI))  ! gpqrs


        totcolsTKI = 0
        do j = 1,nconnectedI
           NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
           NSOMOI = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha+1,NSOMOI+1,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha+1,NSOMOI+1,extype,pmodel,qmodel,2)
           do l = 1,rowsTKI
              do m = 1,colsikpq
                 TKI(l,totcolsTKI+m) = AIJpqContainer(NSOMOalpha+1,NSOMOI+1,extype,pmodel,qmodel,l,m) * psi_coef_config(j)
              enddo
           enddo
           do l = 1,colsikpq
              do m = 1,nconnectedI
                 ! <ij|kl> = (ik|jl)
                 moi = excitationIds(1,j)
                 mok = excitationIds(2,j)
                 moj = excitationIds(1,m)
                 mol = excitationIds(2,m)
                 GIJpqrs(totcolsTKI+l,j) = mo_two_e_integral(moi,moj,mok,mol)
              enddo
           enddo
           totcolsTKI += colsikpq
        end do

        ! Do big BLAS
        ! TODO TKI, size(TKI,1)*size(TKI,2)
        call dgemm('N','N', rowsTKI, nconnectedI, totcolsTKI, 1.d0,  &
          TKI, size(TKI,1), GIJpqrs, size(GIJpqrs,1), 0.d0, &
          TKIGIJ , size(TKIGIJ,1) )

        ! Collect the result
        do j = 1,nconnectedI
           NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
           NSOMOI     = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha+1,NSOMOI+1,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha+1,NSOMOI+1,extype,pmodel,qmodel,2)
           do l = 1,rowsTKI
              do m = 1,colsikpq
                 psi_coef_out(totcolsTKI + m) = AIJpqContainer(NSOMOalpha+1,NSOMOI+1,extype,pmodel,qmodel,l,m) * TKIGIJ(l,j)
              enddo
           enddo
           totcolsTKI += colsikpq
        enddo


     end do
  end do

  end
