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

subroutine calculate_preconditioner_cfg(diag_energies)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Documentation for calculate_preconditioner
  !
  ! Calculates the diagonal energies of
  ! the configurations in psi_configuration
  ! returns : diag_energies :
  END_DOC
  integer :: i,j,k,l,p,q,noccp,noccq, ii, jj
  real*8,intent(out) :: diag_energies(dimBasisCSF)
  integer                            :: nholes
  integer                            :: nvmos
  integer                            :: listvmos(mo_num)
  integer                            :: vmotype(mo_num) ! 1 -> VMO 2 -> SOMO
  integer                            :: listholes(mo_num)
  integer                            :: holetype(mo_num) ! 1-> SOMO 2->DOMO
  integer*8                          :: Idomo
  integer*8                          :: Isomo
  integer*8                          :: Jdomo
  integer*8                          :: Jsomo
  integer*8                          :: diffSOMO
  integer*8                          :: diffDOMO
  integer                            :: NSOMOI
  integer                            :: NSOMOJ
  integer                            :: ndiffSOMO
  integer                            :: ndiffDOMO
  integer                            :: starti, endi, cnti, cntj, rows,cols
  integer                            :: extype,pmodel,qmodel
  integer(bit_kind) :: Icfg(N_INT,2)
  integer(bit_kind) :: Jcfg(N_INT,2)
  integer,external  :: getNSOMO
  real*8, external  :: mo_two_e_integral
  real*8            :: hpp
  real*8            :: meCC
  real*8            :: ecore

  ! initialize energies
  diag_energies = 0.d0

  ! calculate core energy
  call get_core_energy(ecore)
  !diag_energies = ecore

  ! calculate the core energy
  print *,"Core energy=",ref_bitmask_energy

  do i=1,N_configuration

     Isomo = psi_configuration(1,1,i)
     Idomo = psi_configuration(1,2,i)
     Icfg(1,1) = psi_configuration(1,1,i)
     Icfg(1,2) = psi_configuration(1,2,i)
     NSOMOI = getNSOMO(psi_configuration(:,:,i))

     starti = psi_config_data(i,1)
     endi   = psi_config_data(i,2)

     ! find out all pq holes possible
     nholes = 0
     ! holes in SOMO
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
     !do k = 1+n_core_inact_orb,n_core_orb+n_core_inact_act_orb
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
     print *,"I=",i
     call debug_spindet(psi_configuration(1,1,i),N_int)
     call debug_spindet(psi_configuration(1,2,i),N_int)

     do k=1,nholes
        p = listholes(k)
        noccp = holetype(k)

        ! Calculate one-electron
        ! and two-electron coulomb terms
        do l=1,nholes
           q = listholes(l)
           noccq = holetype(l)
           print *,"--------------- K=",p," L=",q

           ! one-electron term
           if(p.EQ.q) then
              hpp = noccq * h_core_ri(p,q)!mo_one_e_integrals(q,q)
           else
              hpp = 0.d0
           endif


           do j=starti,endi
              ! coulomb term
              ! (pp,qq) = <pq|pq>
              if(p.EQ.q) then
                 diag_energies(j) += hpp + 0.5d0 * (noccp * noccq * mo_two_e_integral(p,q,p,q)) !- noccp * mo_two_e_integral(p,q,p,q))
                 print *,"hpp=",hpp,"diga= ",diag_energies(j)
              else
                 diag_energies(j) +=       0.5d0 * noccp * noccq * mo_two_e_integral(p,q,p,q)
                 print *,"diga= ",diag_energies(j)
              endif
           enddo
        enddo


        ! Calculate two-electron
        ! terms type (pk,kp)
        do l=1,nvmos
           !if(.true.)cycle
           q = listvmos(l)
           noccq = vmotype(l)
           if (p.EQ.q) cycle
           print *,"--------------- K=",p," L=",q

           ! ERI term
           ! (pk,kq) = <pq|kk>
           if(noccp .EQ. 1 .AND. noccq .EQ. 0) then
              ! SOMO -> VMO
              NSOMOJ = NSOMOI
              extype = 3
              Jsomo = IBCLR(Isomo,p-1)
              Jsomo = IBSET(Jsomo,q-1)
              Jdomo = Idomo
              Jcfg(1,1) = Jsomo
              Jcfg(1,2) = Jdomo
              pmodel = -1
              qmodel = -1
              call convertOrbIdsToModelSpaceIds(Icfg, Jcfg, p, q, extype, pmodel, qmodel)
              rows = AIJpqMatrixDimsList(NSOMOI,NSOMOJ,extype,pmodel,qmodel,1)
              cols = AIJpqMatrixDimsList(NSOMOI,NSOMOJ,extype,pmodel,qmodel,2)
              call printMatrix(AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,:,:),NBFMax,NBFMax)
              cnti = 1
              cntj = 1
              do ii=1,rows
                 cntj = 1
                 do jj=1,cols
                    meCC = AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,cnti,cntj)
                    meCC *= meCC
                    diag_energies(starti+ii-1) += 0.5d0 * (mo_two_e_integral(p,p,q,q) * meCC)
                    print *,"SOMO->VMO",mo_two_e_integral(p,p,q,q),meCC,NSOMOI,NSOMOJ,"|",rows,cols,">",pmodel,qmodel," diag=",diag_energies(starti+ii-1)
                    cntj += 1
                 enddo
                 cnti += 1
              enddo
           elseif(noccp .EQ. 1 .AND. noccq .EQ. 1) then
              ! SOMO -> SOMO
              NSOMOJ = NSOMOI-2
              extype = 1
              Jsomo = IBCLR(Isomo,p-1)
              Jsomo = IBCLR(Jsomo,q-1)
              Jdomo = IBSET(Idomo,q-1)
              Jcfg(1,1) = Jsomo
              Jcfg(1,2) = Jdomo
              pmodel = -1
              qmodel = -1
              call convertOrbIdsToModelSpaceIds(Icfg, Jcfg, p, q, extype, pmodel, qmodel)
              rows = AIJpqMatrixDimsList(NSOMOI,NSOMOJ,extype,pmodel,qmodel,1)
              cols = AIJpqMatrixDimsList(NSOMOI,NSOMOJ,extype,pmodel,qmodel,2)
              call printMatrix(AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,:,:),NBFMax,NBFMax)
              !print *,Isomo,Idomo,Jsomo,Jdomo,p,q
              !print *,AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,:,:)
              cnti = 1
              cntj = 1
              do ii=1,rows
                 cntj = 1
                 do jj=1,cols
                    meCC = AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,cnti,cntj)
                    meCC *= meCC
                    diag_energies(starti+ii-1) += 0.25d0 * (mo_two_e_integral(p,p,q,q) * meCC )
                    print *,"SOMO->SOMO",mo_two_e_integral(p,p,q,q),meCC,NSOMOI,NSOMOJ,"|",rows,cols,">",pmodel,qmodel," diag=",diag_energies(starti+ii-1)
                    cntj += 1
                 enddo
                 cnti += 1
              enddo
           elseif(noccp .EQ. 2 .AND. noccq .EQ. 0) then
              ! DOMO -> VMO
              NSOMOJ = NSOMOI+2
              extype = 2
              Jsomo = IBSET(Isomo,q-1)
              Jsomo = IBSET(Jsomo,p-1)
              Jdomo = IBCLR(Idomo,p-1)
              Jcfg(1,1) = Jsomo
              Jcfg(1,2) = Jdomo
              pmodel = -1
              qmodel = -1
              call convertOrbIdsToModelSpaceIds(Icfg, Jcfg, p, q, extype, pmodel, qmodel)
              rows = AIJpqMatrixDimsList(NSOMOI,NSOMOJ,extype,pmodel,qmodel,1)
              cols = AIJpqMatrixDimsList(NSOMOI,NSOMOJ,extype,pmodel,qmodel,2)
              call printMatrix(AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,:,:),NBFMax,NBFMax)
              cnti = 1
              cntj = 1
              do ii=1,rows
                 cntj = 1
                 do jj=1,cols
                    meCC = AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,cnti,cntj)
                    meCC *= meCC
                    diag_energies(starti+ii-1) += 0.5d0 * (mo_two_e_integral(p,p,q,q) * meCC)
                    print *,"DOMO->VMO",mo_two_e_integral(p,p,q,q),meCC,NSOMOI,NSOMOJ,"|",rows,cols,">",pmodel,qmodel," diag=",diag_energies(starti+ii-1)
                    cntj += 1
                 enddo
                 cnti += 1
              enddo
           elseif(noccp .EQ. 2 .AND. noccq .EQ. 1) then
              ! DOMO -> SOMO
              NSOMOJ = NSOMOI
              extype = 4
              Jsomo = IBCLR(Isomo,q-1)
              Jsomo = IBSET(Jsomo,p-1)
              Jdomo = IBCLR(Idomo,p-1)
              Jdomo = IBSET(Jdomo,q-1)
              Jcfg(1,1) = Jsomo
              Jcfg(1,2) = Jdomo
              pmodel = -1
              qmodel = -1
              call convertOrbIdsToModelSpaceIds(Icfg, Jcfg, p, q, extype, pmodel, qmodel)
              rows = AIJpqMatrixDimsList(NSOMOI,NSOMOJ,extype,pmodel,qmodel,1)
              cols = AIJpqMatrixDimsList(NSOMOI,NSOMOJ,extype,pmodel,qmodel,2)
              call printMatrix(AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,:,:),NBFMax,NBFMax)
              cnti = 1
              cntj = 1
              do ii=1,rows
                 cntj = 1
                 do jj=1,cols
                    meCC = AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,cnti,cntj)
                    meCC *= meCC
                    diag_energies(starti+ii-1) += 0.5d0 * (mo_two_e_integral(p,p,q,q) * meCC)
                    print *,"DOMO->SOMO",mo_two_e_integral(p,p,q,q),meCC,NSOMOI,NSOMOJ,"|",rows,cols,">",pmodel,qmodel," diag=",diag_energies(starti+ii-1)
                    cntj += 1
                 enddo
                 cnti += 1
              enddo
           else
              print *,"Something is wrong in calculate_preconditioner_cfg"
           endif
        enddo
     enddo
  enddo

end subroutine calculate_preconditioner_cfg

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
  integer(bit_kind) :: alphas_Icfg(N_INT,2,200)
  integer(bit_kind) :: singlesI(N_INT,2,200)
  integer(bit_kind) :: connectedI_alpha(N_INT,2,200)
  integer           :: idxs_singlesI(200)
  integer           :: idxs_connectedI_alpha(200)
  integer(bit_kind) :: psi_configuration_out(N_INT,2,400)
  real*8            :: psi_coef_out(dimBasisCSF)
  logical           :: psi_coef_out_init(dimBasisCSF)
  integer           :: excitationIds_single(2,200)
  integer           :: excitationTypes_single(200)
  integer           :: excitationIds(2,200)
  integer           :: excitationTypes(200)
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
  integer*8 :: MS, Isomo, Idomo, Jsomo, Jdomo
  integer :: moi, moj, mok, mol, starti, endi, startj, endj, cnti, cntj, cntk
  real*8  :: norm_coef_cfg, fac2eints
  real*8  :: norm_coef_det
  real*8  :: meCC1, meCC2
  real*8,dimension(:,:),allocatable :: TKI
  real*8,dimension(:,:),allocatable  :: GIJpqrs
  real*8,dimension(:,:),allocatable  :: TKIGIJ
  real*8, external :: mo_two_e_integral
  real*8, external :: get_two_e_integral
  real*8          :: diag_energies(dimBasisCSF)
  call calculate_preconditioner_cfg(diag_energies)

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
     print *, "i=",i,"coef=",psi_coef_config(i,1), "diagE=",diag_energies(i)
     !call debug_spindet(psi_configuration(1,1,i),N_int)
     !call debug_spindet(psi_configuration(1,2,i),N_int)
     norm_coef_cfg += psi_coef_config(i,1)*psi_coef_config(i,1)
  enddo
  print *,"norm CFG = ",norm_coef_cfg


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
     print *,"-------------------I=",i, nsinglesI
     call debug_spindet(Isomo,N_int)
     call debug_spindet(Idomo,N_int)

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
        call debug_spindet(Jsomo,N_int)
        call debug_spindet(Jdomo,N_int)

        ! Add the hole on J
        if(POPCNT(IAND(Jsomo,IBSET(0_8,q-1))) .EQ. 1) then
           nholes += 1
           listholes(nholes) = q
           holetype(nholes) = 1
        endif
        !if(POPCNT(IAND(Jdomo,IBSET(0_8,q-1))) .EQ. 1 .AND. POPCNT(IAND(Idomo,IBSET(0_8,q-1))).EQ.0) then
        !   nholes += 1
        !   listholes(nholes) = q
        !   holetype(nholes) = 2
        !endif

        print *,"J=",j, "(,",p,q,")", pmodel, qmodel, extype, idxI
        call debug_spindet(singlesI(1,1,j),N_int)
        call debug_spindet(singlesI(1,2,j),N_int)
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
              print *,jj,"sing=",h_core_ri(p,q), meCC1,psi_coef_out(jj)
              cntj += 1
           enddo
           cnti += 1
        enddo

        !!! Two-electron contribution !!!
        do k=1,nholes
           orbk = listholes(k)
           nocck = holetype(k)
           if(k.EQ.p .OR. k.EQ.q) then
              ! 2e integral of the type
              ! (pq,kk) = <pk,qk>
              ! and
              ! (kk,pq) = <kp,kq>
              moi = p
              mok = q
              moj = orbk
              mol = orbk
              !fac2eints *= mo_two_e_integral(moi,moj,mok,mol)
              fac2eints = 0.5d0 * get_two_e_integral(moi,mok,moj,mol,mo_integrals_map) * (nocck)
           else
              ! 2e integral of the type
              ! (pq,kk) = <pk,qk>
              ! and
              ! (kk,pq) = <kp,kq>
              moi = p
              mok = q
              moj = orbk
              mol = orbk
              ! 2 * 0.5 since we do this toice
              ! <I|E_{kk}E_{pq}|J> and <I|E_{pq}E_{kk}|J>
              !fac2eints *= mo_two_e_integral(moi,moj,mok,mol)
              fac2eints = 1.0d0 * get_two_e_integral(moi,mok,moj,mol,mo_integrals_map) * nocck
           endif


           cnti = 1
           cntj = 1
           do ii = starti, endi
              cntj  = 1
              do jj = startj, endj
                 meCC1 = AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,cnti,cntj)
                 psi_coef_out(jj) += meCC1* psi_coef_config(ii,1) * fac2eints
                 psi_coef_out_init(jj) = .True.
                 print *,jj,"doub=",get_two_e_integral(moi,mok,moj,mol,mo_integrals_map),meCC1,psi_coef_out(jj)
                 cntj+=1
              enddo
              cnti += 1
           enddo
           !cnti = 1
           !cntj = 1
           !do ii = starti, endi
           !   cntk  = 1
           !   do kk = startj, endj
           !      meCC1 = AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,cnti,cntk)
           !      cntj = 1
           !      do jj = starti, endi
           !         if(jj.EQ.ii) then
           !            cntj+=1
           !            cycle
           !         endif
           !         meCC2 = AIJpqContainer(NSOMOI,NSOMOJ,extype,pmodel,qmodel,cntj,cntk)
           !         psi_coef_out(jj) += meCC1*meCC2 * psi_coef_config(ii,1) * fac2eints
           !         psi_coef_out_init(jj) = .True.
           !         cntj+=1
           !      enddo
           !      cntk += 1
           !   enddo
           !   cnti += 1
           !enddo
        enddo
        ! Undo setting in listholes
        if(POPCNT(IAND(Jsomo,IBSET(0_8,q-1))) .EQ. 1) then
           nholes -= 1
        endif
        !if(POPCNT(IAND(Jdomo,IBSET(0_8,q-1))) .EQ. 1) then
        !   nholes -= 1
        !endif
     enddo
  enddo

  ! Add the diagonal contribution
  print *,"Done singles"
  do i = 1,dimBasisCSF
     print *, "i=",i,"coef=",psi_coef_config(i,1),psi_coef_out(i)," ini?=",psi_coef_out_init(i)
  enddo

  !!! Double Excitations !!!

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
        print *,"Kalpha=",k
        call debug_spindet(alphas_Icfg(1,1,k),N_int)
        call debug_spindet(alphas_Icfg(1,2,k),N_int)
        ! Now generate all singly excited with respect to a given alpha CFG
        call obtain_connected_I_foralpha(i,alphas_Icfg(:,:,k),connectedI_alpha,idxs_connectedI_alpha,nconnectedI,excitationIds,excitationTypes)

        print *,k,"----> nconnected = ",nconnectedI
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
            print *,"----------------Jcfg------- Isingle=",j
            call debug_spindet(connectedI_alpha(1,1,j),N_int)
            call debug_spindet(connectedI_alpha(1,2,j),N_int)
            print *,"----------------",NSOMOalpha,NSOMOI,"ex=",extype,pmodel,qmodel,"(",rowsikpq,colsikpq,")"
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
           if(p.EQ.q) then
              NSOMOalpha = NSOMOI
              print *,"something is wrong in sigma-vector algorithm.", NSOMOI, NSOMOalpha
           endif
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
                 print *,"integrals (",m,l,")",moi,mok,moj,mol
                 GIJpqrs(totcolsTKI+m,l) = 0.5d0*mo_two_e_integral(moi,moj,mok,mol) ! g(pq,sr) = <ps,qr>
              enddo
           enddo
           totcolsTKI += colsikpq
        enddo


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
           if(p.EQ.q) then
              NSOMOalpha = NSOMOI
              print *,"something is wrong in sigma-vector algorithm.", NSOMOI, NSOMOalpha
           endif
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

     enddo ! loop over alphas
  enddo ! loop over I


  ! Add the diagonal contribution
  do i = 1,dimBasisCSF
     print *, "i=",i,"coef=",psi_coef_config(i,1),psi_coef_out(i)," ini?=",psi_coef_out_init(i)
     !psi_coef_out(i) += 0.0d0*diag_energies(i)*psi_coef_config(i,1)
     psi_coef_out(i) = diag_energies(i)*psi_coef_config(i,1)
     print *, "i=",i,"coef=",psi_coef_out(i)
  enddo

  integer::N_st_loc,startdet,enddet,countdet,ndetI
  real*8 ::psi_energy_loc(1)
  double precision ::psi_s2_loc(N_det,1)
  real*8 ::psi_energy_loc2
  double precision ::psi_coef_out_loc2(N_det,1)
  real*8 :: coefcontrib
  real*8 :: energy_hpsi, energy_qp2, norm_coef_loc
  double precision :: hij
  integer(bit_kind)::tmp_det(N_int)
  integer(bit_kind)::tmp_det2(N_int)
  integer(bit_kind)::tmp_tmp2det(N_int,2)
  integer(bit_kind)::tmp_tmp2det2(N_int,2)
  N_st_loc=1
  psi_energy_loc2=0.d0
  !call u_0_H_u_0(psi_energy_loc2,psi_s2_loc,psi_coef,N_det,psi_det,N_int,N_st_loc,psi_det_size)
  call H_S2_u_0_nstates_openmp(psi_coef_out_loc2,psi_s2_loc,psi_coef,1,N_det)

  psi_coef_out_det = 0.d0

  call convertWFfromCSFtoDET(psi_coef_out,psi_coef_out_det)
  !print *,"energy=",psi_energy_loc2," psi_s2=",psi_s2_loc
  energy_hpsi=0.d0
  energy_qp2=0.d0
  norm_coef_det=0.d0
  norm_coef_loc=0.d0
  countdet=1
  do i = 1,N_configuration
     startdet = psi_configuration_to_psi_det(1,i)
     enddet = psi_configuration_to_psi_det(2,i)
     ndetI = enddet-startdet+1

     do k=1,ndetI
        !norm_coef_det += psi_coef_out_det(countdet,1)*psi_coef_out_det(countdet,1)
        norm_coef_det += psi_coef(countdet,1)*psi_coef(countdet,1)
        norm_coef_loc += psi_coef_out_loc2(countdet,1)*psi_coef_out_loc2(countdet,1)
        energy_qp2 += psi_coef_out_loc2(countdet,1)*psi_coef(countdet,1)
        energy_hpsi += psi_coef_out_det(countdet,1)*psi_coef(countdet,1)
        !print *, "i=",i,ndetI," > ",psi_coef_out_det(startdet+k-1,1)," >> ",psi_coef_out_loc2(startdet+k-1,1)-psi_coef_out_det(startdet+k-1,1)
        print *, "i=",i,ndetI," > ",psi_coef_out_det(startdet+k-1,1)," >> ",psi_coef_out_loc2(startdet+k-1,1)
     enddo
     countdet += ndetI
  enddo
  norm_coef_det = sqrt(norm_coef_det)
  norm_coef_loc = sqrt(norm_coef_loc)
  print *,"norm = ",norm_coef_det, " size=",N_det, " Energy=",energy_hpsi, " Energyqp=",energy_qp2 !+nuclear_repulsion
  print *,"norm = ",norm_coef_det, " size=",N_det, " Energy=",energy_hpsi/norm_coef_det, " Energyqp=",energy_qp2/norm_coef_det !+nuclear_repulsion

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
      real*8          :: diag_energies(dimBasisCSF)
      real*8          :: psi_coef_cfg_out(dimBasisCSF,1)
      real*8          :: psi_coef_det_out(n_det,1)
      integer         :: s, bfIcfg, countcsf
      integer*8         :: Ialpha, Ibeta, Isomo
      !call calculate_preconditioner_cfg(diag_energies)
      !do i=1,N_configuration
      !   print *,i,">",diag_energies(i)
      !enddo
      !call calculate_sigma_vector_cfg(psi_coef_out_det)
      normcfg = 0.d0
      normdet = 0.d0
      call convertWFfromDETtoCSF(psi_coef,psi_coef_cfg_out)
      countcsf = 1
      do i=1,N_configuration
         s = 0
         do k=1,N_int
            if (psi_configuration(k,1,i) == 0_bit_kind) cycle
            s = s + popcnt(psi_configuration(k,1,i))
         enddo
         bfIcfg = max(1,nint((binom(s,(s+1)/2)-binom(s,((s+1)/2)+1))))

         do j = 1,bfIcfg
            print *,countcsf,">",psi_coef_cfg_out(countcsf,1)
            normcfg += psi_coef_cfg_out(countcsf,1)*psi_coef_cfg_out(countcsf,1)
            countcsf += 1
         enddo

      enddo
      call convertWFfromCSFtoDET(psi_coef_cfg_out,psi_coef_det_out)
      do i=1,n_det
         Ialpha = psi_det(1,1,i)
         Ibeta  = psi_det(1,2,i)
         Isomo = IEOR(Ialpha,Ibeta)
         !print *,i,">",psi_coef_det_out(i,1), psi_coef(i,1)
         print *,i,">",psi_coef_det_out(i,1), psi_coef(i,1), abs(psi_coef_det_out(i,1)-psi_coef(i,1)), POPCNT(Isomo)
         normdet += psi_coef_det_out(i,1)*psi_coef_det_out(i,1)
      enddo
      print *,"Norm cfg = ",normcfg," Norm det=",normdet
      end
