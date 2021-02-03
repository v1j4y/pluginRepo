  BEGIN_PROVIDER [ integer, NSOMOMax]
 &BEGIN_PROVIDER [ integer, NCSFMax]
 &BEGIN_PROVIDER [ integer*8, NMO]
 &BEGIN_PROVIDER [ integer, NBFMax]
 &BEGIN_PROVIDER [ integer, dimBasisCSF]
 &BEGIN_PROVIDER [ integer, maxDetDimPerBF]
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
  integer i,j,k,l
  integer startdet,enddet
  integer ncfg,ncfgprev
  integer NSOMO
  integer detDimperBF
  real*8 :: coeff
  maxDetDimPerBF = 0
  detDimperBF = 0
  ncfgprev = 0
  dimBasisCSF = 0
  do i = 2-iand(elec_alpha_num-elec_beta_num,1), elec_num,2
  ncfg = cfg_seniority_index(i) - ncfgprev
  detDimperBF = max(1,int((binom(i,(i+1)/2)-binom(i,((i+1)/2)+1))))
  dimBasisCSF += ncfg * detDimperBF
  if(detDimperBF > maxDetDimPerBF) maxDetDimPerBF = detDimperBF
  ncfgprev = ncfg
  enddo
  END_PROVIDER

  BEGIN_PROVIDER [ real*8, DetToCSFTransformationMatrix, (NSOMOMax,NBFMax,maxDetDimPerBF)]
 &BEGIN_PROVIDER [ real*8, psi_coef_config, (dimBasisCSF)]
  use cfunctions
  implicit none
  BEGIN_DOC
  ! Documentation for DetToCSFTransformationMatrix
  ! Provides the matrix of transformatons for the
  ! conversion between determinant to CSF basis (in BFs)
  END_DOC
  integer*8 :: Isomo
  integer   :: rows, cols, i, j, k
  integer   :: startdet, enddet
  integer*8 MS
  real*8    :: tempBuffer(NBFMax,maxDetDimPerBF)
  real*8    :: tempCoeff(maxDetDimPerBF)
  MS = elec_alpha_num - elec_beta_num
  do i = 2-iand(elec_alpha_num-elec_beta_num,1), elec_num,2
    Isomo = IBSET(0_8, i) - 1_8
    ! rows = Ncsfs
    ! cols = Ndets
    call getCSFtoDETTransformationMatrix(Isomo, MS, NBFMax, maxDetDimPerBF, tempBuffer)
    DetToCSFTransformationMatrix(i,:,:) =  tempBuffer
  enddo

  integer s, bfIcfg
  integer countcsf
  countcsf = 0
  integer countdet
  countdet = 0
  integer istate
  istate = 1
  do i = 1,N_configuration
      startdet = psi_configuration_to_psi_det(1,i)
      enddet = psi_configuration_to_psi_det(2,i)

      do j = startdet, enddet
        tempCoeff(countdet) = psi_coef(psi_configuration_to_psi_det_data(j), istate)
        countdet += 1
      enddo


      s = 0
      do k=1,N_int
        if (psi_configuration(k,1,i) == 0_bit_kind) cycle
        s = s + popcnt(psi_configuration(k,1,i))
      enddo
      bfIcfg = max(1,int((binom(s,(s+1)/2)-binom(s,((s+1)/2)+1))))

      ! perhaps blocking with CFGs of same seniority
      ! can be more efficient
      tempBuffer = DetToCSFTransformationMatrix(s,:,:)

       call dgemm('N','N', NBFMax, 1, maxDetDimPerBF, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, size(tempCoeff,1), 0.d0, psi_coef_config(countcsf), size(psi_coef_config,1))
      !call dgemv('N', NBFMax, maxDetDimPerBF, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, 1, 0.d0, psi_coef_config(countcsf), 1)

      countcsf += bfIcfg
  enddo

  END_PROVIDER

  BEGIN_PROVIDER [ integer, AIJpqMatrixDimsList, (NSOMOMax+1,NSOMOMax+1,4,NSOMOMax,NSOMOMax,2)]
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
  ! Type
  ! 1. SOMO -> VMO
  do i = 0, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i-2,i+2, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,NSOMOMax
           do l = 1,NSOMOMax
              if(k == l) cycle
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
              AIJpqMatrixDimsList(i+1,j+1,1,k,l,1) = rows
              AIJpqMatrixDimsList(i+1,j+1,1,k,l,2) = cols
              AIJpqMatrixDimsList(i+1,j+1,1,l,k,1) = rows
              AIJpqMatrixDimsList(i+1,j+1,1,l,k,2) = cols
              ! j -> i
              AIJpqMatrixDimsList(j+1,i+1,1,k,l,1) = rows
              AIJpqMatrixDimsList(j+1,i+1,1,k,l,2) = cols
              AIJpqMatrixDimsList(j+1,i+1,1,l,k,1) = rows
              AIJpqMatrixDimsList(j+1,i+1,1,l,k,2) = cols
           end do
        end do
     end do
  end do
  ! Type
  ! 2. DOMO -> VMO
  do i = 0, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i-2,i+2, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,NSOMOMax
           do l = 1,NSOMOMax
              if(k == l) cycle
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
              AIJpqMatrixDimsList(i+1,j+1,2,k,l,1) = rows
              AIJpqMatrixDimsList(i+1,j+1,2,k,l,2) = cols
              AIJpqMatrixDimsList(i+1,j+1,2,l,k,1) = rows
              AIJpqMatrixDimsList(i+1,j+1,2,l,k,2) = cols
              ! j -> i
              AIJpqMatrixDimsList(j+1,i+1,2,k,l,1) = rows
              AIJpqMatrixDimsList(j+1,i+1,2,k,l,2) = cols
              AIJpqMatrixDimsList(j+1,i+1,2,l,k,1) = rows
              AIJpqMatrixDimsList(j+1,i+1,2,l,k,2) = cols
           end do
        end do
     end do
  end do
  ! Type
  ! 3. DOMO -> VMO
  do i = 0, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i-2,i+2, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,NSOMOMax
           do l = 1,NSOMOMax
              if(k == l) cycle
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
              AIJpqMatrixDimsList(i+1,j+1,3,k,l,1) = rows
              AIJpqMatrixDimsList(i+1,j+1,3,k,l,2) = cols
              AIJpqMatrixDimsList(i+1,j+1,3,l,k,1) = rows
              AIJpqMatrixDimsList(i+1,j+1,3,l,k,2) = cols
              ! j -> i
              AIJpqMatrixDimsList(j+1,i+1,3,k,l,1) = rows
              AIJpqMatrixDimsList(j+1,i+1,3,k,l,2) = cols
              AIJpqMatrixDimsList(j+1,i+1,3,l,k,1) = rows
              AIJpqMatrixDimsList(j+1,i+1,3,l,k,2) = cols
           end do
        end do
     end do
  end do
  ! Type
  ! 4. DOMO -> SOMO
  do i = 0, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i-2,i+2, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,NSOMOMax
           do l = 1,NSOMOMax
              if(k == l) cycle
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
              AIJpqMatrixDimsList(i+1,j+1,4,k,l,1) = rows
              AIJpqMatrixDimsList(i+1,j+1,4,k,l,2) = cols
              AIJpqMatrixDimsList(i+1,j+1,4,l,k,1) = rows
              AIJpqMatrixDimsList(i+1,j+1,4,l,k,2) = cols
              ! j -> i
              AIJpqMatrixDimsList(j+1,i+1,4,k,l,1) = rows
              AIJpqMatrixDimsList(j+1,i+1,4,k,l,2) = cols
              AIJpqMatrixDimsList(j+1,i+1,4,l,k,1) = rows
              AIJpqMatrixDimsList(j+1,i+1,4,l,k,2) = cols
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
  integer*8 Isomo, Jsomo, tmpsomo
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
  ! Type
  ! 1. SOMO -> SOMO
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i-2,i-2, 2
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        do k = 1,NSOMOMax
           do l = 1,NSOMOMax

              ! Define Jsomo
              Jsomo = IBCLR(Isomo, k-1)
              Jsomo = IBCLR(Jsomo, l-1)

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
  ! Type
  ! 2. DOMO -> VMO
  do i = 2, NSOMOMax, 2
     tmpsomo = ISHFT(1_8,i)-1
     do j = i+2,i+2, 2
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        do k = 1,NSOMOMax
           do l = 1,NSOMOMax
              Isomo = IBCLR(tmpsomo,i-1)
              Isomo = IBCLR(Isomo,j-1)

              ! Define Jsomo
              Jsomo = ISHFT(1_8,i)-1;

              AIJpqContainer(i,j,2,k,l,:,:) = 0.0d0
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
                    AIJpqContainer(i,j,2,k,l,ri,ci) = meMatrix(ri, ci)
                 end do
              end do
           end do
        end do
     end do
  end do
  ! Type
  ! 3. SOMO -> VMO
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i,i, 2
        Jsomo = ISHFT(1_8,i)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        do k = 1,NSOMOMax
           do l = 1,NSOMOMax
              Isomo = ISHFT(1_8,i+1)-1
              Isomo = IBCLR(Isomo,j)
              Jsomo = ISHFT(1_8,i+1)-1
              Jsomo = IBCLR(Jsomo,i)-1

              AIJpqContainer(i,j,3,k,l,:,:) = 0.0d0
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
                    AIJpqContainer(i,j,3,k,l,ri,ci) = meMatrix(ri, ci)
                 end do
              end do
           end do
        end do
     end do
  end do
  ! Type
  ! 4. DOMO -> SOMO
  do i = 2, NSOMOMax, 2
     do j = i,i, 2
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        do k = 1,NSOMOMax
           do l = 1,NSOMOMax
              Isomo = ISHFT(1_8,i+1)-1
              Isomo = IBCLR(Isomo,i)
              Jsomo = ISHFT(1_8,i+1)-1
              Jsomo = IBCLR(Jsomo,j)-1

              AIJpqContainer(i,j,4,k,l,:,:) = 0.0d0
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
                    AIJpqContainer(i,j,4,k,l,ri,ci) = meMatrix(ri, ci)
                 end do
              end do
           end do
        end do
     end do
  end do
  END_PROVIDER
