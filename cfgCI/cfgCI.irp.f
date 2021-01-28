      subroutine printMatrix(mat, rows, cols)
      implicit none
      BEGIN_DOC
      ! Print a 2D matrix
      END_DOC
      integer i,j
      integer*8,intent(in) :: rows
      integer*8,intent(in) :: cols
      real*8,dimension(:,:),intent(in) :: mat(rows,cols)
      print *,""
      do i=1,rows
         do j=1,cols
            write(*,'(F4.2) ',advance="no") mat(i,j)
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
      integer i, j
      integer orbp, orbq
      integer*8 rows
      integer*8 cols
      integer*8 MS
      print *, N_int
      call debug_det(psi_det(1,1,1),N_int)
      call debug_det(psi_det(1,1,2),N_int)
      call debug_spindet(psi_det(1,1,1),N_int)
      call debug_spindet(psi_det(1,1,2),N_int)
      print *,N_int
      cc = "Hello"
      print *,psi_det_size
      do i = 1, 4
         print *, psi_configuration(1,1,i), psi_configuration(1,2,i)
      end do

      integer Nint
      integer(bit_kind), dimension(1,2,100) :: singles
      integer n_singles
      Nint=1
      do i = 1, 1
         call generate_all_singles_cfg(psi_configuration(:,:,i), singles,&
         n_singles, Nint)
         print *,"Number of singles=",n_singles
         do j = 1, 2
            print *, psi_configuration(1,1,i), singles(1,1,j)
            MS = 0
            rows=-1
            cols=-1
            call getApqIJMatrixDims(psi_configuration(1,1,i),           &
                                    singles(1,1,j), &
                                    MS,                       &
                                    rows,                     &
                                    cols)
                                    print *, i,">",rows, cols

         end do
!        call printCFGlist(Nint, n_singles, singles)
      end do

      integer startDet, endDet
      do i = 1, 4
         startDet = psi_configuration_to_psi_det(1,i)
         endDet = psi_configuration_to_psi_det(2,i)
         do j = startDet, endDet
            print *, "\t",i, j, psi_configuration_to_psi_det_data(j)
         end do
      end do
      print *, 'Now starting to read my provider for dims'
      do i = 4,6,2
         do j = i-2,i+2,2
            print *,i,j,AIJpqMatrixDimsList(i,j,1,i,j,1), AIJpqMatrixDimsList(i,j,1,i,j,2)
         end do
      end do
      print *, 'Now starting to read my provider for matrix'
      do i = 4,6,2
         do j = i-2,i+2,2
            rows = AIJpqMatrixDimsList(i,j,1,i,j,1)
            cols = AIJpqMatrixDimsList(i,j,1,i,j,2)
            print *,i,j
            call printMatrix(AIJpqContainer(i,j,1,i,j,:,:),rows,cols)
         end do
      end do
      print *, 'Hello world Tangled with two blocks'
      end
