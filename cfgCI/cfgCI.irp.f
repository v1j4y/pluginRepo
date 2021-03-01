subroutine get_core_energy(ecore)
  implicit none
  BEGIN_DOC
  ! Documentation for get_core_energy
  !
  ! Supplies the core+inactive energy
  END_DOC
  integer           :: i,j,ii,jj
  real*8, external  :: mo_two_e_integral
  real*8,intent(out)      :: ecore
  ecore = 0.d0
  do i=1,n_core_inact_orb
     ii=list_core_inact(i)
     ecore    +=2.D0*mo_one_e_integrals(ii,ii)
     do j=1,n_core_inact_orb
        jj=list_core_inact(j)
        ecore    +=2.D0*mo_two_e_integral(ii,j,ii,j)-mo_two_e_integral(ii,j,jj,i)
    end do
 end do
end subroutine get_core_energy

subroutine calculate_sigma_vector_cfg(psi_coef_out_det)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Documentation for sigma-vector calculation
  !
  ! Calculates the result of the
  ! application of the hamiltonian to the
  ! wavefunction in CFG basis once
  ! TODO : Things prepare outside this routine
  !  1. Touch the providers for
  !     a. ApqIJ containers
  !     b. DET to CSF transformation matrices
  !  2. DET to CSF transcormation
  !  2. CSF to DET back transcormation
  ! returns : psi_coef_out_det :
  END_DOC
  real*8,intent(out):: psi_coef_out_det(N_det,1)
  integer(bit_kind) :: Icfg(N_INT,2)
  integer :: i,j,k,l,p,q,noccp,noccq, ii, jj, m, n, idxI, kk, nocck,orbk
  integer(bit_kind) :: alphas_Icfg(N_INT,2,400)
  integer(bit_kind) :: singlesI(N_INT,2,400)
  integer(bit_kind) :: connectedI_alpha(N_INT,2,400)
  integer           :: idxs_singlesI(400)
  integer           :: idxs_connectedI_alpha(400)
  integer(bit_kind) :: psi_configuration_out(N_INT,2,400)
  real*8            :: psi_coef_out(n_CSF)
  logical           :: psi_coef_out_init(n_CSF)
  integer           :: excitationIds_single(2,400)
  integer           :: excitationTypes_single(400)
  integer           :: excitationIds(2,400)
  integer           :: excitationTypes(400)
  real*8            :: diagfactors(400)
  integer           :: nholes
  integer           :: nvmos
  integer           :: listvmos(mo_num)
  integer           :: vmotype(mo_num) ! 1 -> VMO 2 -> SOMO
  integer           :: listholes(mo_num)
  integer           :: holetype(mo_num) ! 1-> SOMO 2->DOMO
  integer  :: Nalphas_Icfg, nconnectedI, rowsikpq, colsikpq, nsinglesI
  integer  :: extype,NSOMOalpha,NSOMOI,NSOMOJ,pmodel,qmodel
  integer :: getNSOMO
  integer :: totcolsTKI
  integer :: rowsTKI
  integer :: noccpp
  integer*8 :: MS, Isomo, Idomo, Jsomo, Jdomo, Ialpha, Ibeta
  integer :: moi, moj, mok, mol, starti, endi, startj, endj, cnti, cntj, cntk
  real*8  :: norm_coef_cfg, fac2eints
  real*8  :: norm_coef_det
  real*8  :: meCC1, meCC2, diagfac
  real*8,dimension(:,:),allocatable :: TKI
  real*8,dimension(:,:),allocatable  :: GIJpqrs
  real*8,dimension(:,:),allocatable  :: TKIGIJ
  real*8, external :: mo_two_e_integral
  real*8, external :: get_two_e_integral
  real*8          :: diag_energies(n_CSF)
  call calculate_preconditioner_cfg(diag_energies)

  MS = 0
  norm_coef_cfg=0.d0

  psi_coef_out=0.d0
  psi_coef_out_init = .False.

  print *,"CSF basis dim=",n_CSF


  !!! Single Excitations !!!
  do i=1,N_configuration

     Icfg(1,1) = psi_configuration(1,1,i)
     Icfg(1,2) = psi_configuration(1,2,i)
     Isomo = Icfg(1,1)
     Idomo = Icfg(1,2)
     NSOMOI = getNSOMO(Icfg)

     ! find out all pq holes possible
     nholes = 0
     ! holes in SOMO
     ! list_act
     ! list_core
     ! list_core_inact
     ! bitmasks
     !do k = n_core_orb+1,n_core_orb + n_act_orb
     do k = 1,mo_num
        if(POPCNT(IAND(Isomo,IBSET(0_8,k-1))) .EQ. 1) then
           nholes += 1
           listholes(nholes) = k
           holetype(nholes) = 1
        endif
     enddo
     ! holes in DOMO
     !do k = n_core_orb+1,n_core_orb + n_act_orb
     do k = 1,mo_num
        if(POPCNT(IAND(Idomo,IBSET(0_8,k-1))) .EQ. 1) then
           nholes += 1
           listholes(nholes) = k
           holetype(nholes) = 2
        endif
     enddo

     ! find vmos
     listvmos = -1
     vmotype = -1
     nvmos = 0
     !do k = n_core_orb+1,n_core_orb + n_act_orb
     do k = 1,mo_num
        !print *,i,IBSET(0,i-1),POPCNT(IAND(Isomo,(IBSET(0_8,i-1)))), POPCNT(IAND(Idomo,(IBSET(0_8,i-1))))
        if(POPCNT(IAND(Isomo,(IBSET(0_8,k-1)))) .EQ. 0 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,k-1)))) .EQ. 0) then
           nvmos += 1
           listvmos(nvmos) = k
           vmotype(nvmos) = 0
        else if(POPCNT(IAND(Isomo,(IBSET(0_8,k-1)))) .EQ. 1 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,k-1)))) .EQ. 0 ) then
           nvmos += 1
           listvmos(nvmos) = k
           vmotype(nvmos) = 1
        end if
     enddo


     ! Icsf ids
     starti = psi_config_data(i,1)
     endi   = psi_config_data(i,2)
     NSOMOI = getNSOMO(Icfg)
     !print *,"I=",i
     !call debug_spindet(Icfg(1,1),N_int)
     !call debug_spindet(Icfg(1,2),N_int)

     call generate_all_singles_cfg_with_type(Icfg,singlesI,idxs_singlesI,excitationIds_single, &
          excitationTypes_single,nsinglesI,N_int)
     !print *,"-------------------I=",i, nsinglesI, " nholes=",nholes
     !call debug_spindet(Isomo,N_int)
     !call debug_spindet(Idomo,N_int)

     do j = 1,nsinglesI
        idxI = idxs_singlesI(j)
        NSOMOJ = getNSOMO(singlesI(:,:,j))
        p = excitationIds_single(1,j)
        q = excitationIds_single(2,j)
        extype = excitationTypes_single(j)
        ! Off diagonal terms
        call convertOrbIdsToModelSpaceIds(Icfg, singlesI(:,:,j), p, q, extype, pmodel, qmodel)
        Jsomo = singlesI(1,1,j)
        Jdomo = singlesI(1,2,j)
        !call debug_spindet(Jsomo,N_int)
        !call debug_spindet(Jdomo,N_int)

        ! Add the hole on J
        if(POPCNT(IAND(Jsomo,IBSET(0_8,q-1))) .EQ. 1  .AND. POPCNT(IAND(Isomo,IBSET(0_8,q-1))) .EQ. 0) then
           nholes += 1
           listholes(nholes) = q
           holetype(nholes) = 1
        endif
        if((POPCNT(IAND(Jdomo,IBSET(0_8,q-1))) .EQ. 1 .AND. POPCNT(IAND(Idomo,IBSET(0_8,q-1))) .EQ. 0) .AND. POPCNT(IAND(Isomo,IBSET(0_8,q-1))) .EQ. 0) then
           nholes += 1
           listholes(nholes) = q
           holetype(nholes) = 2
        endif

        !print *,"J=",j, "(,",p,q,")", pmodel, qmodel, extype, idxI
        !call debug_spindet(singlesI(1,1,j),N_int)
        !call debug_spindet(singlesI(1,2,j),N_int)
        startj = psi_config_data(idxI,1)
        endj   = psi_config_data(idxI,2)

        !!! One-electron contribution !!!
        cnti = 1
        do ii = starti, endi
           cntj  = 1
           do jj = startj, endj
              meCC1 = AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,cnti,cntj)
              psi_coef_out(jj) += meCC1 * psi_coef_config(ii,1) * h_core_ri(p,q)
              psi_coef_out_init(jj) = .True.
              !print *,jj,"sing=",h_core_ri(p,q), meCC1,psi_coef_config(ii,1),"=",psi_coef_out(jj)
              cntj += 1
           enddo
           cnti += 1
        enddo

        ! Undo setting in listholes
        if(POPCNT(IAND(Jsomo,IBSET(0_8,q-1))) .EQ. 1  .AND. POPCNT(IAND(Isomo,IBSET(0_8,q-1))) .EQ. 0) then
           nholes -= 1
        endif
        if((POPCNT(IAND(Jdomo,IBSET(0_8,q-1))) .EQ. 1 .AND. POPCNT(IAND(Idomo,IBSET(0_8,q-1))) .EQ. 0) .AND. POPCNT(IAND(Isomo,IBSET(0_8,q-1))) .EQ. 0) then
           nholes -= 1
        endif
     enddo
  enddo

  !print *,"Done singles"
  !do i = 1,n_CSF
  !   print *, "i=",i,"coef=",psi_coef_config(i,1),psi_coef_out(i)," ini?=",psi_coef_out_init(i)
  !enddo

  !!! Double Excitations !!!

  ! Loop over all selected configurations
  do i = 1,N_configuration

     Icfg(1,1) = psi_configuration(1,1,i)
     Icfg(1,2) = psi_configuration(1,2,i)

     ! Returns all unique (checking the past) singly excited cfgs connected to I
     call obtain_associated_alphaI(i, Icfg, alphas_Icfg, Nalphas_Icfg)
     ! TODO : remove doubly excited for return
     !print *,i,"Nalphas = ",Nalphas_Icfg
     ! Here we do 2x the loop. One to count for the size of the matrix, then we compute.
     do k = 1,Nalphas_Icfg
        !print *,"Kalpha=",k
        !call debug_spindet(alphas_Icfg(1,1,k),N_int)
        !call debug_spindet(alphas_Icfg(1,2,k),N_int)
        ! Now generate all singly excited with respect to a given alpha CFG
        call obtain_connected_I_foralpha(i,alphas_Icfg(:,:,k),connectedI_alpha,idxs_connectedI_alpha,nconnectedI,excitationIds,excitationTypes,diagfactors)

        !print *,k,"----> nconnected = ",nconnectedI
        if(nconnectedI .EQ. 0) then
           cycle
           !print *,"something is wrong in sigma-vector nconnectedI=",nconnectedI
        endif
        totcolsTKI = 0
        rowsTKI = -1
        do j = 1,nconnectedI
           NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
           NSOMOI = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
           ! for E_pp E_rs and E_ppE_rr case
           if(p.EQ.q) then
              NSOMOalpha = NSOMOI
           endif
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
           !print *,k, Nalphas_Icfg, "idxI=",idxs_connectedI_alpha(j)
           !call debug_spindet(alphas_Icfg(1,1,k),N_int)
           !call debug_spindet(alphas_Icfg(1,2,k),N_int)
           !print *,"----------------Jcfg------- Isingle=",j
           !call debug_spindet(connectedI_alpha(1,1,j),N_int)
           !call debug_spindet(connectedI_alpha(1,2,j),N_int)
           !print *,"----------------",NSOMOalpha,NSOMOI,"ex=",extype,pmodel,qmodel,"(",rowsikpq,colsikpq,")"
        enddo

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

        totcolsTKI = 0
        do j = 1,nconnectedI
           NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
           NSOMOI = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
           ! for E_pp E_rs and E_ppE_rr case
           !call debug_spindet(alphas_Icfg(:,1,k),1)
           !call debug_spindet(alphas_Icfg(:,2,k),1)
           !print *,"det I"
           !call debug_spindet(connectedI_alpha(:,1,j),1)
           !call debug_spindet(connectedI_alpha(:,2,j),1)
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,2)
           !call printMatrix(AIJpqContainer(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,1:rowsikpq,1:colsikpq),rowsikpq,colsikpq)
           !print *,"j=",j,">",rowsikpq,colsikpq,"ex=",extype,"pmod(p)=",p,"qmod(q)=",q," somoI=",NSOMOI," somoa=",NSOMOalpha, " coef=",psi_coef_config(idxs_connectedI_alpha(j),1)
           do l = 1,rowsTKI
              do m = 1,colsikpq
                 TKI(l,totcolsTKI+m) = AIJpqContainer(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,l,m) * psi_coef_config(idxs_connectedI_alpha(j)+m-1,1)
              enddo
           enddo
           do m = 1,colsikpq
              do l = 1,nconnectedI
                 ! <ij|kl> = (ik|jl)
                 moi = excitationIds(1,j) ! p
                 mok = excitationIds(2,j) ! q
                 moj = excitationIds(2,l) ! s
                 mol = excitationIds(1,l) ! r
                 if(moi.EQ.mok .AND. moj.EQ.mol)then
                    diagfac = diagfactors(j)
                    diagfac *= diagfactors(l)
                    !print *,"integrals (",totcolsTKI+m,l,")",mok,moi,mol,moj, "|", diagfac
                    GIJpqrs(totcolsTKI+m,l) = diagfac*0.5d0*mo_two_e_integral(mok,mol,moi,moj) ! g(pq,sr) = <ps,qr>
                 else
                       diagfac = diagfactors(j)*diagfactors(l)
                       !print *,"integrals (",totcolsTKI+m,l,")",mok,moi,mol,moj, "|", diagfac
                       GIJpqrs(totcolsTKI+m,l) = diagfac*0.5d0*mo_two_e_integral(mok,mol,moi,moj) ! g(pq,sr) = <ps,qr>
                    !endif
                 endif
              enddo
           enddo
           totcolsTKI += colsikpq
        enddo


        !print *,"TKI matrix dims= (",rowsTKI,",", totcolsTKI,")"
        !call printMatrix(TKI,rowsTKI,totcolsTKI)
        !print *,"GIJpqrs matrix dims= (",totcolsTKI,",", nconnectedI,")"
        !call printMatrix(GIJpqrs,totcolsTKI,nconnectedI)

        ! Do big BLAS
        ! TODO TKI, size(TKI,1)*size(TKI,2)
        call dgemm('N','N', rowsTKI, nconnectedI, totcolsTKI, 1.d0,  &
             TKI, size(TKI,1), GIJpqrs, size(GIJpqrs,1), 0.d0, &
             TKIGIJ , size(TKIGIJ,1) )

        !print *,"TKIGIJ matrix dims= (",rowsTKI,",", nconnectedI,")"
        !call printMatrix(TKIGIJ,rowsTKI,nconnectedI)

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
                 !print *,"j=",j,pmodel,qmodel,extype," idx=",idxs_connectedI_alpha(j)+m-1, " coef=",AIJpqContainer(NSOMOalpha,NSOMOI,extype,pmodel,qmodel,l,m) * TKIGIJ(l,j), " afc=", AIJpqContainer(Nsomoalpha,NSOMOI,extype,pmodel,qmodel,l,m)," tot=", psi_coef_out(idxs_connectedI_alpha(j)+m-1)
              enddo
           enddo
           totcolsTKI += colsikpq
        enddo

        deallocate(TKI) ! coefficients of CSF
        ! Initialize the inegral container
        ! dims : (totcolsTKI, nconnectedI)
        deallocate(GIJpqrs)  ! gpqrs
        deallocate(TKIGIJ)  ! gpqrs

     enddo ! loop over alphas
  enddo ! loop over I


  ! Add the diagonal contribution
  do i = 1,n_CSF
     !print *, "i=",i,"coef=",psi_coef_config(i,1),psi_coef_out(i)," ini?=",psi_coef_out_init(i)
     psi_coef_out(i) += 1.0d0*diag_energies(i)*psi_coef_config(i,1)
     !psi_coef_out(i) = diag_energies(i)*psi_coef_config(i,1)
     !print *, "i=",i,"coef=",psi_coef_out(i)
  enddo

  integer::N_st_loc,startdet,enddet,countdet,ndetI,ndontmatch
  real*8 ::psi_energy_loc(1)
  double precision ::psi_s2_loc(N_det,1)
  real*8 ::psi_energy_loc2
  double precision ::psi_coef_out_loc2(N_det,1)
  real*8 :: coefcontrib, sqrt2
  real*8 :: energy_hpsi, energy_qp2, norm_coef_loc
  double precision :: hij
  logical :: issame
  integer(bit_kind)::tmp_det(N_int)
  integer(bit_kind)::tmp_det2(N_int)
  integer(bit_kind)::tmp_tmp2det(N_int,2)
  integer(bit_kind)::tmp_tmp2det2(N_int,2)
  N_st_loc=1
  ndontmatch = 0
  psi_energy_loc2=0.d0
  !call u_0_H_u_0(psi_energy_loc2,psi_s2_loc,psi_coef,N_det,psi_det,N_int,N_st_loc,psi_det_size)
  !call H_S2_u_0_nstates_openmp(psi_coef_out_loc2,psi_s2_loc,psi_coef,1,N_det)
  call H_u_0_nstates_openmp(psi_coef_out_loc2,psi_coef,1,N_det)

  psi_coef_out_det = 0.d0

  call convertWFfromCSFtoDET(1,psi_coef_out,psi_coef_out_det)
  ! calculate H|Psi> manually
  !psi_coef_out_det = 0.d0
  !psi_coef_out_det(1,1) = 2.0d0 * h_core_ri(1,1) + 0.0d0 * h_core_ri(2,2)
  !moi = 1
  !mok = 1
  !moj = 2
  !mol = 2
  !psi_coef_out_det(1,1) += 0.5d0 * ( 4.0d0 * mo_two_e_integral(moi,mok,mok,moi) + 0.0d0 * mo_two_e_integral(moi,moj,mok,mol))
  !psi_coef_out_det(1,1) += 0.5d0 * ( 0.0d0 * mo_two_e_integral(moj,moi,mol,mok) + 0.0d0 * mo_two_e_integral(moj,moj,mol,mol))
  !psi_coef_out_det(1,1) += 0.5d0 * ( 2.0d0 * mo_two_e_integral(moi,moj,mol,mok) + 0.0d0 * mo_two_e_integral(moj,moi,mok,mol))
  !psi_coef_out_det(1,1) += 0.5d0 * ( 2.0d0 * mo_two_e_integral(moi,mok,moj,mol) + 0.0d0 * mo_two_e_integral(moj,mol,moi,mok))
  !psi_coef_out_det(2,1) = 0.0d0 * h_core_ri(1,1) + 2.0d0 * h_core_ri(2,2)
  !moi = 1
  !mok = 1
  !moj = 2
  !mol = 2
  !psi_coef_out_det(2,1) += 0.5d0 * ( 0.0d0 * mo_two_e_integral(moi,mok,mok,moi) + 0.0d0 * mo_two_e_integral(moi,moj,mok,mol))
  !psi_coef_out_det(2,1) += 0.5d0 * ( 0.0d0 * mo_two_e_integral(moj,moi,mol,mok) + 4.0d0 * mo_two_e_integral(moj,moj,mol,mol))
  !psi_coef_out_det(2,1) += 0.5d0 * ( 0.0d0 * mo_two_e_integral(moi,moj,mol,mok) + 2.0d0 * mo_two_e_integral(moj,moi,mok,mol))
  !psi_coef_out_det(2,1) += 0.5d0 * ( 0.0d0 * mo_two_e_integral(moi,mok,moj,mol) + 2.0d0 * mo_two_e_integral(moj,mol,moi,mok))
  !print *,"energy=",psi_energy_loc2," psi_s2=",psi_s2_loc
  energy_hpsi=0.d0
  energy_qp2=0.d0
  norm_coef_det=0.d0
  norm_coef_loc=0.d0
  countdet=1
  !print *,"(14,44)=",1.0*mo_two_e_integral(1,4,4,4)*sqrt(2.d0)
  !print *,"(12,24)=",mo_two_e_integral(1,3,4,3)*sqrt(2.d0)
  !print *,"(34,13)=",mo_two_e_integral(3,1,3,4)*sqrt(2.d0)
  sqrt2 = dsqrt(2.0d0)
  !psi_coef_out_det(8,1) =+1.d0*h_core_ri(4,1)
  !psi_coef_out_det(8,1)+= 0.5d0*(-1.d0*mo_two_e_integral(3,4,1,3)*1.d0 &
  !                              +1.d0*mo_two_e_integral(2,4,1,2)*1.d0 &
  !                              +2.d0*mo_two_e_integral(4,2,2,1)*1.d0 &
  !                              +1.d0*(mo_two_e_integral(1,4,1,1)*1.d0 + mo_two_e_integral(2,4,2,1)*1.d0 + mo_two_e_integral(3,4,3,1)*2.d0) &
  !                                    +mo_two_e_integral(2,4,2,1)*1.d0 + mo_two_e_integral(3,4,3,1)*2.d0 + mo_two_e_integral(4,4,4,1)*1.d0)
  !psi_coef_out_det(8,1) = psi_coef_out_det(8,1)/sqrt2
  print *,"1->",h_core_ri(1,4)
  !print *,"1->",-0.5*sqrt2*mo_two_e_integral(1,3,3,4)*1.d0
  !print *,"1->",+0.5*sqrt2*mo_two_e_integral(2,1,4,2)*1.d0
  !print *,"1->",-0.5*sqrt2*mo_two_e_integral(1,1,1,4)*1.d0
  !print *,"1->",-0.5*sqrt2*mo_two_e_integral(2,1,2,4)*2.d0
  !print *,"1->",-0.5*sqrt2*mo_two_e_integral(4,1,4,4)*1.d0
  !print *,"1->",-0.5*sqrt2*mo_two_e_integral(1,2,4,2)*2.d0
  !print *,"1->",-0.5*sqrt2*mo_two_e_integral(1,4,4,4)*2.d0
  do i = 1,N_configuration
     startdet = psi_configuration_to_psi_det(1,i)
     enddet = psi_configuration_to_psi_det(2,i)
     ndetI = enddet-startdet+1

     do k=1,ndetI
        Ialpha= psi_det(1,1,startdet+k-1)
        Ibeta = psi_det(1,2,startdet+k-1)
        Isomo = IEOR(Ialpha,Ibeta)
        Idomo = IAND(Ialpha,Ibeta)
        !norm_coef_det += psi_coef_out_det(countdet,1)*psi_coef_out_det(countdet,1)
        norm_coef_det += psi_coef(countdet,1)*psi_coef(countdet,1)
        norm_coef_loc += psi_coef_out_loc2(countdet,1)*psi_coef_out_loc2(countdet,1)
        energy_qp2 += psi_coef_out_loc2(countdet,1)*psi_coef(countdet,1)
        energy_hpsi += psi_coef_out_det(countdet,1)*psi_coef(countdet,1)
        issame = .False.
        !if(abs(abs(psi_coef_out_loc2(startdet+k-1,1))-abs(psi_coef_out_det(startdet+k-1,1))) .LT. 1.0e-8) issame = .True.
        if(abs(psi_coef_out_loc2(startdet+k-1,1)-psi_coef_out_det(startdet+k-1,1)) .LT. 1.0e-5) then
           issame = .True.
           print *, "i=",i,countdet,POPCNT(Isomo), startdet+k-1," > ",psi_coef_out_det(startdet+k-1,1)," >> ",psi_coef_out_loc2(startdet+k-1,1)," |", issame
        else
           call debug_spindet(Isomo,1)
           call debug_spindet(Idomo,1)
           print *, "i=",i,countdet,POPCNT(Isomo), startdet+k-1," > ",psi_coef_out_det(startdet+k-1,1)," >> ",psi_coef_out_loc2(startdet+k-1,1)," |", issame
           ndontmatch +=1
        endif
        !print *, "i=",i,ndetI," > ",psi_coef_out_det(startdet+k-1,1)," >> ",psi_coef_out_loc2(startdet+k-1,1)
     enddo
     countdet += ndetI
  enddo
  norm_coef_det = sqrt(norm_coef_det)
  norm_coef_loc = sqrt(norm_coef_loc)
  print *,"dont match=",ndontmatch,"norm = ",norm_coef_det, " size=",N_det, " Energy=",energy_hpsi/norm_coef_det, " Energyqp=",energy_qp2/norm_coef_det !+nuclear_repulsion

end subroutine calculate_sigma_vector

subroutine calculate_sigma_vector_cfg_test(psi_coef_out_det)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Documentation for sigma-vector calculation
  !
  ! Calculates the result of the
  ! application of the hamiltonian to the
  ! wavefunction in CFG basis once
  ! TODO : Things prepare outside this routine
  !  1. Touch the providers for
  !     a. ApqIJ containers
  !     b. DET to CSF transformation matrices
  !  2. DET to CSF transcormation
  !  2. CSF to DET back transcormation
  ! returns : psi_coef_out_det :
  END_DOC
  real*8,intent(out):: psi_coef_out_det(N_det,1)
  integer(bit_kind) :: Icfg(N_INT,2)
  integer :: i,j,k,l,p,q,noccp,noccq, ii, jj, m, n, idxI, kk, nocck,orbk
  integer(bit_kind) :: alphas_Icfg(N_INT,2,400)
  integer(bit_kind) :: singlesI(N_INT,2,400)
  integer(bit_kind) :: connectedI_alpha(N_INT,2,400)
  integer           :: idxs_singlesI(400)
  integer           :: idxs_connectedI_alpha(400)
  integer(bit_kind) :: psi_configuration_out(N_INT,2,400)
  real*8            :: psi_coef_out(n_CSF,1)
  logical           :: psi_coef_out_init(n_CSF)
  integer           :: excitationIds_single(2,400)
  integer           :: excitationTypes_single(400)
  integer           :: excitationIds(2,400)
  integer           :: excitationTypes(400)
  real*8            :: diagfactors(400)
  integer           :: nholes
  integer           :: nvmos
  integer           :: listvmos(mo_num)
  integer           :: vmotype(mo_num) ! 1 -> VMO 2 -> SOMO
  integer           :: listholes(mo_num)
  integer           :: holetype(mo_num) ! 1-> SOMO 2->DOMO
  integer  :: Nalphas_Icfg, nconnectedI, rowsikpq, colsikpq, nsinglesI
  integer  :: extype,NSOMOalpha,NSOMOI,NSOMOJ,pmodel,qmodel
  integer :: getNSOMO
  integer :: totcolsTKI
  integer :: rowsTKI
  integer :: noccpp
  integer*8 :: MS, Isomo, Idomo, Jsomo, Jdomo, Ialpha, Ibeta
  integer :: moi, moj, mok, mol, starti, endi, startj, endj, cnti, cntj, cntk
  real*8  :: norm_coef_cfg, fac2eints
  real*8  :: norm_coef_det
  real*8  :: meCC1, meCC2, diagfac
  real*8,dimension(:,:),allocatable :: TKI
  real*8,dimension(:,:),allocatable  :: GIJpqrs
  real*8,dimension(:,:),allocatable  :: TKIGIJ
  real*8, external :: mo_two_e_integral
  real*8, external :: get_two_e_integral
  real*8          :: diag_energies(n_CSF)

  !touch dettocsftransformationmatrix psi_coef_config psi_config_data psi_csf_to_config_data


  MS = 0
  norm_coef_cfg=0.d0

  psi_coef_out=0.d0
  psi_coef_out_init = .False.

  print *,"CSF basis dim=",n_CSF

    
  call calculate_sigma_vector_cfg_nst(psi_coef_out, psi_coef_config, 1, n_CSF, 1, n_CSF, 0, 1)

  ! Add the diagonal contribution
  !do i = 1,n_CSF
  !   !print *, "i=",i,"coef=",psi_coef_config(i,1),psi_coef_out(i)," ini?=",psi_coef_out_init(i)
  !   psi_coef_out(i) += 1.0d0*diag_energies(i)*psi_coef_config(i,1)
  !   !psi_coef_out(i) = diag_energies(i)*psi_coef_config(i,1)
  !   !print *, "i=",i,"coef=",psi_coef_out(i)
  !enddo

  integer::N_st_loc,startdet,enddet,countdet,ndetI,ndontmatch
  real*8 ::psi_energy_loc(1)
  double precision ::psi_s2_loc(N_det,1)
  real*8 ::psi_energy_loc2
  double precision ::psi_coef_out_loc2(N_det,1)
  real*8 :: coefcontrib, sqrt2
  real*8 :: energy_hpsi, energy_qp2, norm_coef_loc
  double precision :: hij
  logical :: issame
  integer(bit_kind)::tmp_det(N_int)
  integer(bit_kind)::tmp_det2(N_int)
  integer(bit_kind)::tmp_tmp2det(N_int,2)
  integer(bit_kind)::tmp_tmp2det2(N_int,2)
  N_st_loc=1
  ndontmatch = 0
  energy_qp2=0.d0
  psi_energy_loc2=0.d0
  !call u_0_H_u_0(psi_energy_loc2,psi_s2_loc,psi_coef,N_det,psi_det,N_int,N_st_loc,psi_det_size)
  call H_u_0_nstates_openmp(psi_coef_out_loc2,psi_coef,1,N_det)
  do i=1,N_det
      energy_qp2 += psi_coef_out_loc2(i,1)*psi_coef(i,1)
  enddo
 double precision :: i_H_psi_array(N_states)

 energy_qp2=0.d0
 norm_coef_loc = 0.d0
 do i=1,N_det
  call i_H_psi(psi_det(1,1,i), psi_det, psi_coef, N_int, N_det, &
               size(psi_coef,1), N_states, i_H_psi_array)
  do j=1,1
    norm_coef_loc += psi_coef(i,j)*psi_coef(i,j)
    energy_qp2 += i_H_psi_array(j) * psi_coef(i,j)
  enddo
 enddo

 print *, 'Energy:'
 do i=1,1
   print *, energy_qp2/norm_coef_loc
 enddo

  psi_coef_out_det = 0.d0

  call convertWFfromCSFtoDET(1,psi_coef_out(:,1),psi_coef_out_det)
  energy_hpsi=0.d0
  norm_coef_det=0.d0
  norm_coef_loc=0.d0
  countdet=1
  sqrt2 = dsqrt(2.0d0)
  print *,"1->",h_core_ri(1,4)
  do i = 1,N_configuration
     startdet = psi_configuration_to_psi_det(1,i)
     enddet = psi_configuration_to_psi_det(2,i)
     ndetI = enddet-startdet+1

     do k=1,ndetI
        Ialpha= psi_det(1,1,startdet+k-1)
        Ibeta = psi_det(1,2,startdet+k-1)
        Isomo = IEOR(Ialpha,Ibeta)
        Idomo = IAND(Ialpha,Ibeta)
        !norm_coef_det += psi_coef_out_det(countdet,1)*psi_coef_out_det(countdet,1)
        norm_coef_det += psi_coef(startdet+k-1,1)*psi_coef(startdet+k-1,1)
        norm_coef_loc += psi_coef_out_loc2(startdet+k-1,1)*psi_coef_out_loc2(startdet+k-1,1)
        !energy_qp2 += psi_coef_out_loc2(startdet+k-1,1)*psi_coef(startdet+k-1,1)
        energy_hpsi += psi_coef_out_det(startdet+k-1,1)*psi_coef(startdet+k-1,1)
        issame = .False.
        !if(abs(abs(psi_coef_out_loc2(startdet+k-1,1))-abs(psi_coef_out_det(startdet+k-1,1))) .LT. 1.0e-8) issame = .True.
        if(abs(psi_coef_out_loc2(startdet+k-1,1)-psi_coef_out_det(startdet+k-1,1)) .LT. 1.0e-5) then
           issame = .True.
           print *, "i=",i,countdet,POPCNT(Isomo), startdet+k-1," > ",psi_coef_out_det(startdet+k-1,1)," >> ",psi_coef_out_loc2(startdet+k-1,1)," |", issame
        else
           call debug_spindet(Isomo,1)
           call debug_spindet(Idomo,1)
           print *, "i=",i,countdet,POPCNT(Isomo), startdet+k-1," > ",psi_coef_out_det(startdet+k-1,1)," >> ",psi_coef_out_loc2(startdet+k-1,1)," |", issame
           ndontmatch +=1
        endif
        !print *, "i=",i,ndetI," > ",psi_coef_out_det(startdet+k-1,1)," >> ",psi_coef_out_loc2(startdet+k-1,1)
     countdet += 1
     enddo
  enddo
  norm_coef_det = sqrt(norm_coef_det)
  norm_coef_loc = sqrt(norm_coef_loc)
  print *,"dont match=",ndontmatch,"norm = ",norm_coef_det, " size=",N_det, " Energy=",energy_hpsi +nuclear_repulsion, " Energyqp=",energy_qp2 +nuclear_repulsion, "Nuclear=",nuclear_repulsion

end subroutine calculate_sigma_vector

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
            write(*,'(F20.8,4X)',advance="no") mat(i,j)
         enddo
         print *,""
      enddo
      end subroutine printMatrix

      program cfgCI
      use cfunctions
      implicit none
      BEGIN_DOC
!     TODO : Put the documentation of the program here
      END_DOC
      integer         :: i,j,k,l,p,q
      real*8          :: normcfg, normdet
      real*8          :: psi_coef_out_det(N_det,1)
      real*8          :: diag_energies(n_CSF)
      real*8          :: psi_coef_cfg_out(n_CSF,1)
      real*8          :: psi_coef_det_out(n_det,1)
      integer         :: s, bfIcfg, countcsf
      integer*8         :: Ialpha, Ibeta, Isomo
      !call calculate_preconditioner_cfg(diag_energies)
      !do i=1,N_configuration
      !   print *,i,">",diag_energies(i)
      !enddo
      !call calculate_sigma_vector_cfg(psi_coef_out_det)
      call calculate_sigma_vector_cfg_test(psi_coef_out_det)
      ! Testing CSF->DET->CSF
      !normcfg = 0.d0
      !normdet = 0.d0
      !call convertWFfromDETtoCSF(psi_coef,psi_coef_cfg_out)
      !countcsf = 1
      !do i=1,N_configuration
      !   s = 0
      !   do k=1,N_int
      !      if (psi_configuration(k,1,i) == 0_bit_kind) cycle
      !      s = s + popcnt(psi_configuration(k,1,i))
      !   enddo
      !   bfIcfg = max(1,nint((binom(s,(s+1)/2)-binom(s,((s+1)/2)+1))))

      !   do j = 1,bfIcfg
      !      print *,countcsf,">",psi_coef_cfg_out(countcsf,1)
      !      normcfg += psi_coef_cfg_out(countcsf,1)*psi_coef_cfg_out(countcsf,1)
      !      countcsf += 1
      !   enddo

      !enddo
      !call convertWFfromCSFtoDET(psi_coef_cfg_out,psi_coef_det_out)
      !do i=1,n_det
      !   Ialpha = psi_det(1,1,i)
      !   Ibeta  = psi_det(1,2,i)
      !   Isomo = IEOR(Ialpha,Ibeta)
      !   !print *,i,">",psi_coef_det_out(i,1), psi_coef(i,1)
      !   print *,i,">",psi_coef_det_out(i,1), psi_coef(i,1), abs(psi_coef_det_out(i,1)-psi_coef(i,1)), POPCNT(Isomo)
      !   normdet += psi_coef_det_out(i,1)*psi_coef_det_out(i,1)
      !enddo
      !print *,"Norm cfg = ",normcfg," Norm det=",normdet
      end
