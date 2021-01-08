program cfgCI
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  character*32 cc
  print *, N_int
  print *, psi_configuration(1,1,:200)
  print *, psi_configuration(1,2,:200)
  call debug_det(psi_det(1,1,1),N_int)
  call debug_det(psi_det(1,1,2),N_int)
  call debug_spindet(psi_det(1,1,1),N_int)
  call debug_spindet(psi_det(1,1,2),N_int)
  print *,N_int
  cc = "Hello"
  print *,psi_det_size
  call printCFGlist(cfg_seniority_index,N_int,psi_det_size)
  print *, 'Hello world'
end
