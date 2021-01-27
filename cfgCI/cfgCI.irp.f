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
      do i = 1, 15
         print *, psi_configuration(1,1,i), psi_configuration(1,2,i)
      end do
      call printCFGlist(N_int, psi_det_size, psi_configuration)

      integer Nint
      integer(bit_kind), dimension(1,2,100) :: singles
      integer n_singles
      Nint=1
      do i = 1, 1
         call generate_all_singles_cfg(psi_configuration(:,:,i), singles,&
         n_singles, Nint)
         print *,"Number of singles=",n_singles
         do j = 1, 20
            print *, psi_configuration(1,2,i), singles(1,2,j)
            MS = 0
            rows=-1
            cols=-1
            call getApqIJMatrixDims(singles(1,2,j),           &
                                    psi_configuration(1,2,i), &
                                    MS,                       &
                                    rows,                     &
                                    cols)
                                    print *, i,">",rows, cols

         end do
         call printCFGlist(Nint, n_singles, singles)
      end do
      print *, 'Hello world Tangled'
      end
