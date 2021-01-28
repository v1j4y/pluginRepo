module cfunctions
        use, intrinsic :: ISO_C_BINDING
      interface
         subroutine printcfglist(nint, ncfgs, cfglist) bind(C, name='printCFGList')
           import C_INT32_T, C_INT64_T
         integer(kind=C_INT32_T) :: nint
         integer(kind=C_INT32_T) :: ncfgs
         integer(kind=C_INT64_T) :: cfglist(nint,2,ncfgs)
       end subroutine printcfglist
      end interface
      interface
         subroutine getApqIJMatrixDims(Isomo, Jsomo, MS, rowsout, colsout) &
              bind(C, name='getApqIJMatrixDims')
           import C_INT32_T, C_INT64_T
           integer(kind=C_INT64_T),value,intent(in) :: Isomo ! CSFI
           integer(kind=C_INT64_T),value,intent(in) :: Jsomo ! CSFJ
           integer(kind=C_INT64_T),value,intent(in) :: MS    ! Ms = 2*Spin
           integer(kind=C_INT64_T),intent(out):: rowsout
           integer(kind=C_INT64_T),intent(out):: colsout
         end subroutine getApqIJMatrixDims
      end interface
      interface
         subroutine getApqIJMatrixDriver(Isomo, Jsomo, orbp, orbq,  &
              MS, NMO, CSFICSFJApqIJ, rowsout, colsout) bind(C, name='getApqIJMatrixDriverArrayInp')
           import C_INT32_T, C_INT64_T, C_DOUBLE
           integer(kind=C_INT64_T),value,intent(in) :: Isomo
           integer(kind=C_INT64_T),value,intent(in) :: Jsomo
           integer(kind=C_INT32_T),value,intent(in) :: orbp
           integer(kind=C_INT32_T),value,intent(in) :: orbq
           integer(kind=C_INT64_T),value,intent(in) :: MS
           integer(kind=C_INT64_T),value,intent(in) :: NMO
           real(kind=C_DOUBLE),intent(out) :: CSFICSFJApqIJ(200,200)
           integer(kind=C_INT64_T),intent(out) :: rowsout
           integer(kind=C_INT64_T),intent(out) :: colsout
           !integer(kind=C_INT32_T),dimension(rowApqIJ,colApqIJ) :: ApqIJ
         end subroutine getApqIJMatrixDriver
      end interface
    end module cfunctions
