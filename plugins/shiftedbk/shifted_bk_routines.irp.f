use selection_types

 BEGIN_PROVIDER [ double precision, fock_diag_tmp_, (2,mo_tot_num+1,Nproc) ]
&BEGIN_PROVIDER [ integer, current_generator_, (Nproc) ]
&BEGIN_PROVIDER [ double precision, a_h_i, (N_det, Nproc) ]
&BEGIN_PROVIDER [ double precision, a_s2_i, (N_det, Nproc) ]
&BEGIN_PROVIDER [ type(selection_buffer), sb, (Nproc) ]
&BEGIN_PROVIDER [ double precision, N_det_increase_factor ]
  implicit none
  integer :: i
  integer :: n_det_add

  N_det_increase_factor = 1d0

  current_generator_(:) = 0
  do i=1,Nproc
    n_det_add = max(1, int(float(N_det) * N_det_increase_factor))
    call create_selection_buffer(n_det_add, n_det_add*2, sb(i))
  end do
  a_h_i = 0d0
  a_s2_i = 0d0
 END_PROVIDER



subroutine delta_ij_done()
  implicit none
  integer :: i, n_det_add
  
  call sort_selection_buffer(sb(1))

  do i=2,Nproc
    call sort_selection_buffer(sb(i))
    call merge_selection_buffers(sb(i), sb(1))
  end do
  
  call sort_selection_buffer(sb(1))

  call fill_H_apply_buffer_no_selection(sb(1)%cur,sb(1)%det,N_int,0) 
  call copy_H_apply_buffer_to_wf()
  if (s2_eig.or.(N_states > 1) ) then
    call make_s2_eigenfunction
  endif
  !call save_wavefunction
  n_det_add = max(1, int(float(N_det) * N_det_increase_factor))
  do i=1,Nproc
    call delete_selection_buffer(sb(i))
    call create_selection_buffer(n_det_add, n_det_add*2, sb(i))
  end do
  !delta_ij = 0d0
end subroutine


subroutine dress_with_alpha_buffer(Nstates,Ndet,Nint,delta_ij_loc, i_gen, minilist, det_minilist, n_minilist, alpha, iproc)
 use bitmasks
 implicit none
  BEGIN_DOC
  !delta_ij_loc(:,:,1) : dressing column for H
  !delta_ij_loc(:,:,2) : dressing column for S2
  !i_gen : generator index in psi_det_generators
  !minilist : indices of determinants connected to alpha ( in psi_det_sorted )
  !n_minilist : size of minilist
  !alpha : alpha determinant
  END_DOC
  integer, intent(in)             :: Nint, Ndet, Nstates, n_minilist, iproc, i_gen
  integer(bit_kind), intent(in)   :: alpha(Nint,2), det_minilist(Nint, 2, n_minilist)
  integer,intent(in)              :: minilist(n_minilist)
  double precision, intent(inout) :: delta_ij_loc(Nstates,N_det,2)
  double precision :: haa, hij, sij
  double precision, external :: diag_H_mat_elem_fock
  integer                        :: i,j,k,l,m, l_sd
  double precision :: hdress, sdress
  double precision :: de, a_h_psi(Nstates), c_alpha, contrib


  a_h_psi = 0d0
  
  if(current_generator_(iproc) /= i_gen) then
    current_generator_(iproc) = i_gen
    call build_fock_tmp(fock_diag_tmp_(1,1,iproc),psi_det_generators(1,1,i_gen),N_int)
  end if

  haa = diag_H_mat_elem_fock(psi_det_generators(1,1,i_gen),alpha,fock_diag_tmp_(1,1,iproc),N_int)

  do l_sd=1,n_minilist
    call i_h_j_s2(alpha,det_minilist(1,1,l_sd),N_int,hij, sij)
    a_h_i(l_sd, iproc) = hij
    a_s2_i(l_sd, iproc) = sij
    do i=1,Nstates
      a_h_psi(i) += hij * psi_coef(minilist(l_sd), i)
    end do
  end do

  contrib = 0d0

  do i=1,Nstates
    de = E0_denominator(i) - haa
    if(DABS(de) < 1D-5) cycle

    c_alpha = a_h_psi(i) / de
    contrib = min(contrib, c_alpha * a_h_psi(i))
    
    do l_sd=1,n_minilist
      hdress = c_alpha * a_h_i(l_sd, iproc)
      sdress = c_alpha * a_s2_i(l_sd, iproc)
      delta_ij_loc(i, minilist(l_sd), 1) += hdress
      delta_ij_loc(i, minilist(l_sd), 2) += sdress
    end do
  end do

  call add_to_selection_buffer(sb(iproc), alpha, contrib)

end subroutine


BEGIN_PROVIDER [ logical, initialize_E0_denominator ]
    implicit none
    BEGIN_DOC
    ! If true, initialize pt2_E0_denominator
    END_DOC
    initialize_E0_denominator = .True.
END_PROVIDER


BEGIN_PROVIDER [ double precision, E0_denominator, (N_states) ]
  implicit none
  BEGIN_DOC
  ! E0 in the denominator of the PT2
  END_DOC
  if (initialize_E0_denominator) then
   E0_denominator(1:N_states) = psi_energy(1:N_states)
 ! call ezfio_get_full_ci_zmq_energy(pt2_E0_denominator(1))
 ! pt2_E0_denominator(1) -= nuclear_repulsion
 ! pt2_E0_denominator(1:N_states) = HF_energy - nuclear_repulsion
 ! pt2_E0_denominator(1:N_states) = barycentric_electronic_energy(1:N_states)
  else
    E0_denominator = -huge(1.d0)
  endif
END_PROVIDER


