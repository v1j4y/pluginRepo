  BEGIN_PROVIDER [ integer, NSOMOMax]
 &BEGIN_PROVIDER [ integer, NCSFMax]
 &BEGIN_PROVIDER [ integer*8, NMO]
 &BEGIN_PROVIDER [ integer, NBFMax]
  implicit none
  BEGIN_DOC
  ! Documentation for NSOMOMax
  ! The maximum number of SOMOs for the current calculation.
  ! required for the calculation of prototype arrays.
  END_DOC
  NSOMOMax = 8
  NBFMax = 14
  NCSFMax = 14 ! TODO: NCSFs for MS=0
  NMO = NSOMOMax ! TODO: remove this
  END_PROVIDER

  BEGIN_PROVIDER [ integer, AIJpqMatrixDimsList, (NSOMOMax,NSOMOMax,4,NSOMOMax,NSOMOMax,2)]
 &BEGIN_PROVIDER [ integer, rowsmax]
 &BEGIN_PROVIDER [ integer, colsmax]
  use cfunctions
  implicit none
  BEGIN_DOC
  ! Documentation for AIJpqMatrixList
  ! The prototype matrix containing the <I|E_{pq}|J>
  ! matrices for each I,J somo pair and orb ids.
  END_DOC
  integer i,j,k,l
  integer*8 Isomo, Jsomo
  Isomo = 0
  Jsomo = 0
  integer rows, cols
  rows = -1
  cols = -1
  integer*8 MS
  MS = 0
  rowsmax = 0
  colsmax = 0
  print *,"NSOMOMax = ",NSOMOMax
  !allocate(AIJpqMatrixDimsList(NSOMOMax,NSOMOMax,4,NSOMOMax,NSOMOMax,2))
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1,i)-1
     do j = i-2,i+2, 2
        Jsomo = ISHFT(1,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,NSOMOMax
           do l = k,NSOMOMax
              call getApqIJMatrixDims(Isomo,           &
                   Jsomo, &
                   MS,                       &
                   rows,                     &
                   cols)
              print *, i,j,k,l,">",Isomo,Jsomo,">",rows, cols
              if(rowsmax .LT. rows) then
                 rowsmax = rows
              end if
              if(colsmax .LT. cols) then
                 colsmax = cols
              end if
              ! i -> j
              AIJpqMatrixDimsList(i,j,1,k,l,1) = rows
              AIJpqMatrixDimsList(i,j,1,k,l,2) = cols
              AIJpqMatrixDimsList(i,j,1,l,k,1) = rows
              AIJpqMatrixDimsList(i,j,1,l,k,2) = cols
              ! j -> i
              AIJpqMatrixDimsList(j,i,1,k,l,1) = rows
              AIJpqMatrixDimsList(j,i,1,k,l,2) = cols
              AIJpqMatrixDimsList(j,i,1,l,k,1) = rows
              AIJpqMatrixDimsList(j,i,1,l,k,2) = cols
           end do
        end do
     end do
  end do
  print *,"Rowsmax=",rowsmax," Colsmax=",colsmax
  END_PROVIDER

  BEGIN_PROVIDER [ real*8, AIJpqContainer, (NSOMOMax,NSOMOMax,4,NSOMOMax+4,NSOMOMax+4,NBFMax,NBFMax)]
  use cfunctions
  implicit none
  BEGIN_DOC
  ! Documentation for AIJpqMatrixList
  ! The prototype matrix containing the <I|E_{pq}|J>
  ! matrices for each I,J somo pair and orb ids.
  !
  ! Due to the different types of excitations which
  ! include DOMOs and VMOs two prototype DOMOs and two
  ! prototype VMOs are needed. Therefore
  ! the 4th and 5th dimensions are NSOMOMax+4 and NSOMOMax+4
  ! respectively.
  !
  ! The type of excitations are ordered as follows:
  ! Type 1 - SOMO -> SOMO
  ! Type 2 - DOMO -> VMO
  ! Type 3 - SOMO -> VMO
  ! Type 4 - DOMO -> SOMO
  END_DOC
  integer i,j,k,l, orbp, orbq, ri, ci
  orbp = 0
  orbq = 0
  integer*8 Isomo, Jsomo
  Isomo = 0
  Jsomo = 0
  integer rows, cols
  rows = -1
  cols = -1
  integer*8 MS
  MS = 0
  real*8,dimension(:,:),allocatable :: meMatrix
  ! allocate matrix
  allocate(meMatrix(rowsmax,colsmax))
  print *,"NSOMOMax = ",NSOMOMax
  !allocate(AIJpqMatrixDimsList(NSOMOMax,NSOMOMax,4,NSOMOMax,NSOMOMax,2))
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(ISHFT(1,i)-1,1)
     do j = i-2,i+2, 2
        Jsomo = ISHFT(1,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        do k = i-1,i+1
           do l = 1,1
              AIJpqContainer(i,j,1,k,l,:,:) = 0.0d0
              call getApqIJMatrixDims(Isomo,           &
                   Jsomo, &
                   MS,                       &
                   rows,                     &
                   cols)

              orbp = k
              orbq = l
              ! fill matrix
              call getApqIJMatrixDriver(Isomo,           &
                   Jsomo, &
                   orbp,                     &
                   orbq,                     &
                   MS,                       &
                   NMO,                      &
                   meMatrix,                 &
                   rows,                     &
                   cols)
             print *, i,j,k,l,">",Isomo,Jsomo,">",rows, cols,">",rowsmax,colsmax
              ! i -> j
             do ri = 1,rows
                 do ci = 1,cols
                    AIJpqContainer(i,j,1,k,l,ri,ci) = meMatrix(ri, ci)
                 end do
              end do
           end do
        end do
     end do
  end do
  END_PROVIDER
