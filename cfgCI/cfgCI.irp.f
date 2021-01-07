program cfgCI
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  character*32 cc
  print *, cfg_seniority_index
  call debug_det(psi_det(1,1,1),N_int)
  call debug_det(psi_det(1,1,2),N_int)
  call debug_spindet(psi_det(1,1,1),N_int)
  call debug_spindet(psi_det(1,1,2),N_int)
  print *,N_int
  cc = "Hello"
  call printCFGlist(cfg_seniority_index,15)
  print *, 'Hello world'
end
