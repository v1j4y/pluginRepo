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
  ! Note that here we need NSOMOMax + 2 sizes
  NBFMax = 42
  NCSFMax = 42 ! TODO: NCSFs for MS=0
  maxDetDimPerBF =  252
  NMO = NSOMOMax ! TODO: remove this
  integer i,j,k,l
  integer startdet,enddet
  integer ncfg,ncfgprev
  integer NSOMO
  integer dimcsfpercfg
  integer detDimperBF
  real*8 :: coeff
  integer MS
  integer ncfgpersomo
  detDimperBF = 0
  MS = elec_alpha_num-elec_beta_num
  ! number of cfgs = number of dets for 0 somos
  dimBasisCSF = cfg_seniority_index(0)-1
  ncfgprev = cfg_seniority_index(0)-1
  do i = 0-iand(MS,1), NSOMOMax-2,2
     if(cfg_seniority_index(i+2) == -1)then
        ncfgpersomo = N_configuration
     else
        ncfgpersomo = cfg_seniority_index(i+2)
     endif
  ncfg = ncfgpersomo - 1 - ncfgprev
  !detDimperBF = max(1,nint((binom(i,(i+1)/2))))
  dimcsfpercfg = max(1,nint((binom(i,(i+1)/2)-binom(i,((i+1)/2)+1))))
  dimBasisCSF += ncfg * dimcsfpercfg
  print *,i,">",ncfg,",",detDimperBF,">",dimcsfpercfg
  !if(detDimperBF > maxDetDimPerBF) maxDetDimPerBF = detDimperBF
  ncfgprev = ncfg
  enddo
  END_PROVIDER

  BEGIN_PROVIDER [ real*8, DetToCSFTransformationMatrix, (0:NSOMOMax,NBFMax,maxDetDimPerBF)]
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
  print *,"Maxbfdim=",NBFMax
  print *,"Maxdetdim=",maxDetDimPerBF
  print *,"dimBasisCSF=",dimBasisCSF
  DetToCSFTransformationMatrix(0,:,:) = 1.0
  do i = 2-iand(elec_alpha_num-elec_beta_num,1), NSOMOMax,2
    Isomo = IBSET(0_8, i) - 1_8
    ! rows = Ncsfs
    ! cols = Ndets
    call getCSFtoDETTransformationMatrix(Isomo, MS, NBFMax, maxDetDimPerBF, tempBuffer)
    DetToCSFTransformationMatrix(i,:,:) =  tempBuffer
  enddo

  integer s, bfIcfg
  integer countcsf
  countcsf = 1
  integer countdet
  countdet = 1
  integer istate
  istate = 1
  do i = 1,N_configuration
      startdet = psi_configuration_to_psi_det(1,i)
      enddet = psi_configuration_to_psi_det(2,i)

      do j = startdet, enddet
        tempCoeff(countdet) = psi_coef(psi_configuration_to_psi_det_data(j), istate)
        countdet += 1
      enddo


      s = 1
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

  BEGIN_PROVIDER [ integer, AIJpqMatrixDimsList, (0:NSOMOMax,0:NSOMOMax,4,NSOMOMax,NSOMOMax,2)]
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
  integer*8 Isomo, Jsomo, tmpsomo
  Isomo = 0
  Jsomo = 0
  integer rows, cols
  rows = -1
  cols = -1
  integer*8 MS
  MS = elec_alpha_num-elec_beta_num
  integer nsomomin
  nsomomin = elec_alpha_num-elec_beta_num
  rowsmax = 0
  colsmax = 0
  print *,"NSOMOMax = ",NSOMOMax
  !allocate(AIJpqMatrixDimsList(NSOMOMax,NSOMOMax,4,NSOMOMax,NSOMOMax,2))
  ! Type
  ! 1. SOMO -> SOMO
  print *,"Doing SOMO->SOMO"
  do i = 2-iand(nsomomin,1), NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i-2,i-2, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LT. 0) then
           cycle
        end if
        do k = 1,i
           do l = 1,i
              ! Define Jsomo
              if(k.NE.l)then
                 Jsomo = IBCLR(Isomo, k-1)
                 Jsomo = IBCLR(Jsomo, l-1)
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif

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
           end do
        end do
     end do
  end do
  ! Type
  ! 2. DOMO -> VMO
  print *,"Doing DOMO->VMO"
  do i = 0+iand(nsomomin,1), NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     tmpsomo = ISHFT(1_8,i+2)-1
     do j = i+2,i+2, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,j
           do l = 1,j
              if(k .NE. l) then
              Isomo = IBCLR(tmpsomo,k-1)
              Isomo = IBCLR(Isomo,l-1)

              ! Define Jsomo
              Jsomo = ISHFT(1_8,i)-1;
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif

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
              AIJpqMatrixDimsList(i,j,2,k,l,1) = rows
              AIJpqMatrixDimsList(i,j,2,k,l,2) = cols
           end do
        end do
     end do
  end do
  ! Type
  ! 3. DOMO -> SOMO
  print *,"Doing DOMO->SOMO"
  do i = 2-iand(nsomomin,1), NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i,i, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,i
           do l = 1,i
              if(k.NE.l)then
              Isomo = ISHFT(1_8,i+1)-1
              Isomo = IBCLR(Isomo,k)
              Jsomo = ISHFT(1_8,j+1)-1
              Jsomo = IBCLR(Jsomo,l)
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif
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
              AIJpqMatrixDimsList(i,j,3,k,l,1) = rows
              AIJpqMatrixDimsList(i,j,3,k,l,2) = cols
           end do
        end do
     end do
  end do
  ! Type
  ! 4. SOMO -> VMO
  print *,"Doing SOMO->VMO"
  do i = 2-iand(nsomomin,1), NSOMOMax, 2
     do j = i,i, 2
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,i
           do l = 1,i
              if(k.NE.l)then
              Isomo = ISHFT(1_8,i+1)-1
              Isomo = IBCLR(Isomo,k)
              Jsomo = ISHFT(1_8,j+1)-1
              Jsomo = IBCLR(Jsomo,l)
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif
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
              AIJpqMatrixDimsList(i,j,4,k,l,1) = rows
              AIJpqMatrixDimsList(i,j,4,k,l,2) = cols
           end do
        end do
     end do
  end do
  print *,"Rowsmax=",rowsmax," Colsmax=",colsmax
  END_PROVIDER

  BEGIN_PROVIDER [ real*8, AIJpqContainer, (0:NSOMOMax,0:NSOMOMax,4,NSOMOMax,NSOMOMax,NBFMax,NBFMax)]
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
  touch AIJpqMatrixDimsList
  real*8,dimension(:,:),allocatable :: meMatrix
  integer maxdim
  maxdim = max(rowsmax,colsmax)
  ! allocate matrix
  allocate(meMatrix(maxdim,maxdim))
  print *,"rowsmax =",rowsmax," colsmax=",colsmax
  print *,"NSOMOMax = ",NSOMOMax
  !allocate(AIJpqMatrixDimsList(NSOMOMax,NSOMOMax,4,NSOMOMax,NSOMOMax,2))
  ! Type
  ! 1. SOMO -> SOMO
  print *,"Doing SOMO -> SOMO"
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i-2,i-2, 2
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        print *,"i,j=",i,j
        do k = 1,i
           do l = 1,i

              ! Define Jsomo
              if(k .NE. l) then
                 Jsomo = IBCLR(Isomo, k-1)
                 Jsomo = IBCLR(Jsomo, l-1)
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif

              print *,"k,l=",k,l
              call debug_spindet(Jsomo,1)
              call debug_spindet(Isomo,1)

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
  print *,"Doing DOMO -> VMO"
  do i = 0, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     tmpsomo = ISHFT(1_8,i+2)-1
     do j = i+2,i+2, 2
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        Jsomo = ISHFT(1_8,j)-1
        do k = 1,j
           do l = 1,j
              if(k .NE. l) then
                 Isomo = IBCLR(tmpsomo,k-1)
                 Isomo = IBCLR(Isomo,l-1)
                 ! Define Jsomo
                 Jsomo = ISHFT(1_8,j)-1;
              else
                 Isomo = ISHFT(1_8,j)-1
                 Isomo = IBCLR(Isomo,1-1)
                 Isomo = IBCLR(Isomo,j-1)
                 Jsomo = ISHFT(1_8,j)-1
              endif

              print *,"k,l=",k,l
              call debug_spindet(Jsomo,1)
              call debug_spindet(Isomo,1)

              !AIJpqContainer(i,j,2,k,l,:,:) = 0.0d0
              call getApqIJMatrixDims(Isomo,           &
                   Jsomo, &
                   MS,                       &
                   rows,                     &
                   cols)

              print *,"Done Dims"
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
              print *,"Done allocate"
           end do
        end do
     end do
  end do
  ! Type
  ! 3. SOMO -> VMO
  print *,"Doing SOMO -> VMO"
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i,i, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        do k = 1,i
           do l = 1,i
              if(k .NE. l) then
              Isomo = ISHFT(1_8,i+1)-1
              Isomo = IBCLR(Isomo,k)
              Jsomo = ISHFT(1_8,j+1)-1
              Jsomo = IBCLR(Jsomo,l)
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif

              print *,"k,l=",k,l
              call debug_spindet(Jsomo,1)
              call debug_spindet(Isomo,1)

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
  print *,"Doing DOMO -> SOMO"
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i,i, 2
        Jsomo = ISHFT(1_8,i)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        do k = 1,i
           do l = 1,i
              if(k .NE. l) then
              Isomo = ISHFT(1_8,i+1)-1
              Isomo = IBCLR(Isomo,k)
              Jsomo = ISHFT(1_8,j+1)-1
              Jsomo = IBCLR(Jsomo,l)
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif

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
