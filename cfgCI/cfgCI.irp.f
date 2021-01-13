program cfgCI
  use cfunctions
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  character*32 cc
  integer i
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
  print *, 'Hello world'
end
