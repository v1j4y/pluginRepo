subroutine obtain_connected_I_foralpha(idxI, Ialpha, connectedI, idxs_connectedI, nconnectedI, excitationIds, excitationTypes)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Documentation for obtain_connected_I_foralpha
  ! This function returns all those selected configurations
  ! which are connected to the input configuration
  ! Ialpha by a single excitation.
  !
  ! The type of excitations are ordered as follows:
  ! Type 1 - SOMO -> SOMO
  ! Type 2 - DOMO -> VMO
  ! Type 3 - SOMO -> VMO
  ! Type 4 - DOMO -> SOMO
  !
  ! Order of operators
  ! \alpha> = a^\dag_p a_q |I> = E_pq |I>
  END_DOC
  integer          ,intent(in)             :: idxI
  integer(bit_kind),intent(in)             :: Ialpha(N_int,2)
  integer(bit_kind),intent(out)            :: connectedI(N_int,2,*)
  integer          ,intent(out)            :: idxs_connectedI(*)
  integer,intent(out)                      :: nconnectedI
  integer,intent(out)                      :: excitationIds(2,*)
  integer,intent(out)                      :: excitationTypes(*)
  integer*8                                :: Idomo
  integer*8                                :: Isomo
  integer*8                                :: Jdomo
  integer*8                                :: Jsomo
  integer*8                                :: IJsomo
  integer*8                                :: diffSOMO
  integer*8                                :: diffDOMO
  integer                                  :: ndiffSOMO
  integer                                  :: ndiffDOMO
  integer :: i,j,k,l,p,q,nsomoJ,nsomoalpha,starti,endi,extyp,nholes
  integer :: listholes(mo_num)
  integer :: holetype(mo_num)

  ! find out all pq holes possible
  nholes = 0
  ! holes in SOMO
  Isomo = psi_configuration(1,1,idxI)
  Idomo = psi_configuration(1,2,idxI)
  do i = n_core_orb+1,n_core_orb + n_act_orb
     if(POPCNT(IAND(Isomo,IBSET(0_8,i-1))) .EQ. 1) then
        nholes += 1
        listholes(nholes) = i
        holetype(nholes) = 1
     endif
  end do
  ! holes in DOMO
  do i = n_core_orb+1,n_core_orb + n_act_orb
     if(POPCNT(IAND(Idomo,IBSET(0_8,i-1))) .EQ. 1) then
        nholes += 1
        listholes(nholes) = i
        holetype(nholes) = 2
     endif
  end do

  nconnectedI = 0

  p = 0
  q = 0
  do i=idxI+1,N_configuration
     Isomo = Ialpha(1,1)
     Idomo = Ialpha(1,2)
     Jsomo = psi_configuration(1,1,i)
     Jdomo = psi_configuration(1,2,i)
     !call debug_spindet(Isomo,1)
     !call debug_spindet(Idomo,1)
     !print *,"-J--i=",i,Idomo,Jdomo,">",N_configuration
     !call debug_spindet(Jsomo,1)
     !call debug_spindet(Jdomo,1)
     diffSOMO = IEOR(Isomo,Jsomo)
     diffDOMO = IEOR(Idomo,Jdomo)
     ndiffSOMO = POPCNT(diffSOMO)
     ndiffDOMO = POPCNT(diffDOMO)
     if((ndiffSOMO + ndiffDOMO) .EQ. 0) cycle
     !print *,"-I--i=",i,diffSOMO,diffDOMO!Isomo,Jsomo,ndiffSOMO,ndiffDOMO
     !print *,POPCNT(IEOR(diffSOMO,diffDOMO)), ndiffDOMO
     if(POPCNT(IEOR(diffSOMO,diffDOMO)) .LE. 1 .AND. ndiffDOMO .LT. 3) then
     !call debug_spindet(Isomo,1)
     !call debug_spindet(Idomo,1)
     !print *,"-J--i=",i,Idomo,Jdomo,">",N_configuration
     !call debug_spindet(Jsomo,1)
     !call debug_spindet(Jdomo,1)
        select case(ndiffDOMO)
        case (0)
           ! SOMO -> VMO
           !print *,"obt SOMO -> VMO"
           extyp = 3
           IJsomo = IEOR(Isomo, Jsomo)
           p = TRAILZ(IAND(Isomo,IJsomo)) + 1
           IJsomo = IBCLR(IJsomo,p-1)
           q = TRAILZ(IJsomo) + 1
        case (1)
           ! DOMO -> VMO
           ! or
           ! SOMO -> SOMO
           nsomoJ = POPCNT(Jsomo)
           nsomoalpha = POPCNT(Isomo)
           if(nsomoJ .GT. nsomoalpha) then
              ! DOMO -> VMO
              !print *,"obt DOMO -> VMO"
              extyp = 2
              p = TRAILZ(IEOR(Idomo,Jdomo)) + 1
              Isomo = IEOR(Isomo, Jsomo)
              Isomo = IBCLR(Isomo,p-1)
              q = TRAILZ(Isomo) + 1
           else
              ! SOMO -> SOMO
              !print *,"obt SOMO -> SOMO"
              extyp = 1
              q = TRAILZ(IEOR(Idomo,Jdomo)) + 1
              Isomo = IEOR(Isomo, Jsomo)
              Isomo = IBCLR(Isomo,q-1)
              p = TRAILZ(Isomo) + 1
           end if
        case (2)
           ! DOMO -> SOMO
           !print *,"obt DOMO -> SOMO"
           extyp = 4
           IJsomo = IEOR(Isomo, Jsomo)
           p = TRAILZ(IAND(Jsomo,IJsomo)) + 1
           IJsomo = IBCLR(IJsomo,p-1)
           q = TRAILZ(IJsomo) + 1
        case default
           print *,"something went wront in get connectedI"
        end select
        starti = psi_config_data(i,1)
        endi   = psi_config_data(i,2)
        nconnectedI += 1
        connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
        idxs_connectedI(nconnectedI)=starti
        excitationIds(1,nconnectedI)=p
        excitationIds(2,nconnectedI)=q
        excitationTypes(nconnectedI) = extyp
        print *,"------ > output p,q in obt=",p,q
     endif
  end do

end subroutine obtain_connected_I_foralpha
