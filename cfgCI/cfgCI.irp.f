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
  integer           :: excitationIds(200,2)
  integer           :: excitationTypes(200)
  integer  :: Nalphas_Icfg, nconnectedI
  MS = 0

  ! CFG: 2 2 2 0 0 0
  ! CFG: 2 2 1 1 0 0
  ! CFG: 2 1 1 1 1 0
  !Icfg(1,1,1) = 0;
  !Icfg(1,2,1) = ISHFT(1,3)-1;
  !Icfg(1,1,2) = ISHFT(3,2);
  !Icfg(1,2,2) = ISHFT(1,2)-1;
  !Icfg(1,1,3) = ISHFT(ISHFT(1,4),1);
  !Icfg(1,2,3) = 1
  !integer Nint;
  !Nint = 1
  !integer n_singles_Icfg
  !! Generate all singles excitations
  !print *,"--------Gen----------"
  !call debug_spindet(Icfg(:,:,1),Nint)
  !call generate_all_singles_cfg_with_type(Icfg(:,:,1),singles_Icfg, &
  !     ex_type_singles, n_singles_Icfg, Nint)
  !print *,"nsingles=",n_singles_Icfg
  !do i = 1,n_singles_Icfg
  !   call debug_spindet(singles_Icfg(1,1,i),Nint)
  !   call debug_spindet(singles_Icfg(1,2,i),Nint)
  !   call getApqIJMatrixDims(Icfg(1,1,1),           &
  !                           singles_Icfg(1,1,i), &
  !                           MS,                       &
  !                           rows,                     &
  !                           cols)
  !   print *,i," ",singles_Icfg(1,1,i), "-",singles_Icfg(1,2,i), ex_type_singles(i), rows,cols
  !end do

  ! Loop over singles
  do i = 1,N_configuration
     Icfg(1,1) = psi_configuration(1,1,i)
     Icfg(1,2) = psi_configuration(1,2,i)
     call obtain_associated_alphaI(i, Icfg, alphas_Icfg, Nalphas_Icfg)
     print *,i,"Nalphas = ",Nalphas_Icfg
     do k = 1,Nalphas_Icfg
        print *,k, Nalphas_Icfg
        call debug_spindet(alphas_Icfg(1,1,k),N_int)
        call debug_spindet(alphas_Icfg(1,2,k),N_int)
        ! Now generate all singly excited with respect to a given alpha CFG
        call obtain_connected_I_foralpha(alphas_Icfg(:,:,k),connectedI_alpha,nconnectedI,excitationIds,excitationTypes)
        print *,k,"----> nconnected = ",nconnectedI
        do j = 1,nconnectedI
        print *,"----------------",i,j,k,nconnectedI
        call debug_spindet(connectedI_alpha(1,1,j),N_int)
        call debug_spindet(connectedI_alpha(1,2,j),N_int)
        end do
     end do
  end do

  end
