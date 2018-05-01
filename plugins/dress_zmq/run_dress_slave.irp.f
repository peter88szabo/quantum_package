use bitmasks

BEGIN_PROVIDER [ integer, fragment_count ]
  implicit none
  BEGIN_DOC
    ! Number of fragments for the deterministic part
  END_DOC
  fragment_count = 1
END_PROVIDER


subroutine run_dress_slave(thread,iproce,energy)
  use f77_zmq
  implicit none

  double precision, intent(in)    :: energy(N_states_diag)
  integer,  intent(in)            :: thread, iproce
  integer                        :: rc, i, subset, i_generator

  integer                        :: worker_id, task_id, ctask, ltask
  character*(5120)                :: task

  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_push_socket
  integer(ZMQ_PTR)               :: zmq_socket_push

  logical :: done

  double precision,allocatable :: dress_detail(:)
  integer :: ind
  
  double precision,allocatable :: delta_ij_loc(:,:,:)
  double precision :: div(N_states) 
  integer :: h,p,n,i_state
  logical :: ok

  integer, allocatable :: int_buf(:)
  double precision, allocatable :: double_buf(:)
  integer(bit_kind), allocatable :: det_buf(:,:,:)
  integer :: N_buf(3)
  logical :: last
  integer, external :: omp_get_thread_num 
  double precision, allocatable :: delta_det(:,:,:,:), cp(:,:,:,:)
  integer :: toothMwen
  logical :: fracted
  double precision :: fac
      


  if(iproce /= 0) stop "RUN DRESS SLAVE is OMP"
  
  allocate(delta_det(N_states, N_det, 0:comb_teeth+1, 2))
  allocate(cp(N_states, N_det, N_cp, 2))
  delta_det = 0d9
  cp = 0d0

  
  task(:) = CHAR(0)

  
  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  zmq_socket_push      = new_zmq_push_socket(thread)
  call connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread)
  if(worker_id == -1) then
    print *, "WORKER -1"
    call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
    call end_zmq_push_socket(zmq_socket_push,thread)
    return
  end if
  do i=1,N_states
    div(i) = psi_coef(dressed_column_idx(i), i)
  end do
  
  integer :: iproc, cur_cp, done_for(0:N_cp)
  integer, allocatable :: tasks(:)
  integer :: lastCp(Nproc)
  integer :: lastSent, lastSendable
  logical :: send
  lastCp = 0
  lastSent = 0
  send = .false.
  done_for = 0

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(int_buf, double_buf, det_buf, delta_ij_loc, task, task_id) &
  !$OMP PRIVATE(toothMwen, fracted, fac) &
  !$OMP PRIVATE(send, i_generator, subset, iproc, N_buf)
  iproc = omp_get_thread_num()+1
  allocate(int_buf(N_dress_int_buffer)) 
  allocate(double_buf(N_dress_double_buffer)) 
  allocate(det_buf(N_int, 2, N_dress_det_buffer)) 
  allocate(delta_ij_loc(N_states,N_det,2)) 
  do
    !$OMP CRITICAL (SENDAGE)
    call get_task_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id, task)
    !$OMP END CRITICAL (SENDAGE)
    task = task//"  0"
    if(task_id == 0) then
      print *, "DONEDONE"
      exit   !! LAST MESSAGE ???
    end if
    read (task,*) subset, i_generator
    
    
    if(done_cp_at_det(i_generator) < lastCp(iproc)) stop 'loop = .false.'
    !$OMP CRITICAL
    send = .false.
    lastSendable = N_cp*2
    do i=1,Nproc
      lastSendable = min(lastCp(iproc), lastSendable)
    end do
    lastSendable -= 1
    if(lastSendable > lastSent) then
      lastSent = lastSendable
      send = .true.
    end if
    !$OMP END CRITICAL
     
    if(send) then
      !$OMP CRITICAL
      N_buf = (/0,1,0/)
  
      delta_ij_loc = 0d0
      cur_cp = lastSent
      if(cur_cp < 1) stop "cur_cp < 1"
      do i=1,cur_cp
        delta_ij_loc(:,:,:) += cp(:,:,i,:)
      end do

      delta_ij_loc(:,:,:) = delta_ij_loc(:,:,:) / cps_N(cur_cp)
      do i=cp_first_tooth(cur_cp)-1,0,-1
       delta_ij_loc(:,:,:) = delta_ij_loc(:,:,:) +delta_det(:,:,i,:)
      end do
      !$OMP END CRITICAL
      !$OMP CRITICAL (SENDAGE)
      call push_dress_results(zmq_socket_push, done_for(cur_cp), cur_cp, delta_ij_loc, int_buf, double_buf, det_buf, N_buf, -1)
      !$OMP END CRITICAL (SENDAGE)
    end if


    !$OMP ATOMIC
    done_for(done_cp_at_det(i_generator)) += 1

    delta_ij_loc(:,:,:) = 0d0
    call generator_start(i_generator, iproc)
    call alpha_callback(delta_ij_loc, i_generator, subset, iproc)
    call generator_done(i_generator, int_buf, double_buf, det_buf, N_buf, iproc)
     
    !if(.false.) then
    !$OMP CRITICAL
    do i=1,N_cp
      fac = cps(i_generator, i) * dress_weight_inv(i_generator) * comb_step
      if(fac == 0d0) cycle
      cp(:,:,i,1) += (delta_ij_loc(:,:,1) * fac)
      cp(:,:,i,2) += (delta_ij_loc(:,:,2) * fac)
    end do


    toothMwen = tooth_of_det(i_generator)
    fracted = (toothMwen /= 0)
    if(fracted) fracted = (i_generator == first_det_of_teeth(toothMwen))
    if(fracted) then
      delta_det(:,:,toothMwen-1, 1) += delta_ij_loc(:,:,1) * (1d0-fractage(toothMwen))
      delta_det(:,:,toothMwen-1, 2) += delta_ij_loc(:,:,2) * (1d0-fractage(toothMwen))
      delta_det(:,:,toothMwen  , 1) += delta_ij_loc(:,:,1) * (fractage(toothMwen))
      delta_det(:,:,toothMwen  , 2) += delta_ij_loc(:,:,2) * (fractage(toothMwen))
    else
      delta_det(:,:,toothMwen  , 1) += delta_ij_loc(:,:,1)
      delta_det(:,:,toothMwen  , 2) += delta_ij_loc(:,:,2)
    end if

      
    !$OMP END CRITICAL
    !end if

    !$OMP CRITICAL (SENDAGE)
    call push_dress_results(zmq_socket_push, i_generator, -1, delta_ij_loc, int_buf, double_buf, det_buf, N_buf, task_id)
    call task_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id)
    !$OMP END CRITICAL (SENDAGE)
    lastCp(iproc) = done_cp_at_det(i_generator) 
  end do
  !$OMP END PARALLEL
  
  call sleep(10) 
  call disconnect_from_taskserver(zmq_to_qp_run_socket,zmq_socket_push,worker_id)
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_push_socket(zmq_socket_push,thread)
end subroutine



subroutine push_dress_results(zmq_socket_push, ind, cur_cp, delta_loc, int_buf, double_buf, det_buf, N_bufi, task_id)
  use f77_zmq
  implicit none

  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  double precision, intent(inout)   :: delta_loc(N_states, N_det, 2)
  real(kind=4), allocatable   :: delta_loc4(:,:,:)
  double precision, intent(in) :: double_buf(*)
  integer, intent(in) :: int_buf(*)
  integer(bit_kind), intent(in) :: det_buf(N_int, 2, *)
  integer, intent(in) :: N_bufi(3)
  integer :: N_buf(3)
  integer, intent(in) :: ind, cur_cp, task_id
  integer :: rc, i, j, k, l
  double precision :: contrib(N_states)
  

  
  rc = f77_zmq_send( zmq_socket_push, ind, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop "push"
  
  rc = f77_zmq_send( zmq_socket_push, cur_cp, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop "push"
  

  if(cur_cp /= -1) then
    rc = f77_zmq_send( zmq_socket_push, delta_loc(1,1,1), 8*N_states*N_det, ZMQ_SNDMORE)
    if(rc /= 8*N_states*N_det) stop "push"

    rc = f77_zmq_send( zmq_socket_push, delta_loc(1,1,2), 8*N_states*N_det, ZMQ_SNDMORE)
    if(rc /= 8*N_states*N_det) stop "push"
  else
    contrib = 0d0
    
    do i=1,N_det
      contrib(:) += delta_loc(:,i, 1) * psi_coef(i, :)
    end do

    rc = f77_zmq_send( zmq_socket_push, contrib, 8*N_states, ZMQ_SNDMORE)
    if(rc /= 8*N_states) stop "push"
    
    N_buf = N_bufi
    N_buf = (/0,1,0/)

    rc = f77_zmq_send( zmq_socket_push, N_buf, 4*3, ZMQ_SNDMORE)
    if(rc /= 4*3) stop "push5" 
    
    if(N_buf(1) > N_dress_int_buffer) stop "run_dress_slave N_buf bad size?"
    if(N_buf(2) > N_dress_double_buffer) stop "run_dress_slave N_buf bad size?"
    if(N_buf(3) > N_dress_det_buffer) stop "run_dress_slave N_buf bad size?"

    
    if(N_buf(1) > 0) then
      rc = f77_zmq_send( zmq_socket_push, int_buf, 4*N_buf(1), ZMQ_SNDMORE)
      if(rc /= 4*N_buf(1)) stop "push6"
    end if
    
    if(N_buf(2) > 0) then
      rc = f77_zmq_send( zmq_socket_push, double_buf, 8*N_buf(2), ZMQ_SNDMORE)
      if(rc /= 8*N_buf(2)) stop "push8"
    end if

    if(N_buf(3) > 0) then
      rc = f77_zmq_send( zmq_socket_push, det_buf, 2*N_int*bit_kind*N_buf(3), ZMQ_SNDMORE)
      if(rc /= 2*N_int*bit_kind*N_buf(3)) stop "push10"
    end if

    rc = f77_zmq_send( zmq_socket_push, task_id, 4, 0)
    if(rc /= 4) stop "push11"
  end if

! Activate is zmq_socket_push is a REQ
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(2) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
IRP_ENDIF

end subroutine


subroutine pull_dress_results(zmq_socket_pull, ind, cur_cp, delta_loc, int_buf, double_buf, det_buf, N_buf, task_id, contrib)
  use f77_zmq
  implicit none
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  integer, intent(out) :: cur_cp
  double precision, intent(inout) :: delta_loc(N_states, N_det, 2)
  double precision, intent(out) :: double_buf(*), contrib(N_states)
  integer, intent(out) :: int_buf(*)
  integer(bit_kind), intent(out) :: det_buf(N_int, 2, *)
  integer, intent(out) :: ind
  integer, intent(out) :: task_id
  integer :: rc, i, j, k
  integer, intent(out) :: N_buf(3)


  
  rc = f77_zmq_recv( zmq_socket_pull, ind, 4, 0)
  if(rc /= 4) stop "pulla"
 
  rc = f77_zmq_recv( zmq_socket_pull, cur_cp, 4, 0)
  if(rc /= 4) stop "pulla"
  
  
  

  if(cur_cp /= -1) then
    rc = f77_zmq_recv( zmq_socket_pull, delta_loc(1,1,1), N_states*8*N_det, 0)
    if(rc /= 8*N_states*N_det) stop "pullc"
  
    rc = f77_zmq_recv( zmq_socket_pull, delta_loc(1,1,2), N_states*8*N_det, 0)
    if(rc /= 8*N_states*N_det) stop "pulld"
  else
    rc = f77_zmq_recv( zmq_socket_pull, contrib, 8*N_states, 0)
    if(rc /= 8*N_states) stop "pullc"
  
    rc = f77_zmq_recv( zmq_socket_pull, N_buf, 4*3, 0)
    if(rc /= 4*3) stop "pull" 
    if(N_buf(1) > N_dress_int_buffer) stop "run_dress_slave N_buf bad size?"
    if(N_buf(2) > N_dress_double_buffer) stop "run_dress_slave N_buf bad size?"
    if(N_buf(3) > N_dress_det_buffer) stop "run_dress_slave N_buf bad size?"

    
    if(N_buf(1) > 0) then
      rc = f77_zmq_recv( zmq_socket_pull, int_buf, 4*N_buf(1), 0)
      if(rc /= 4*N_buf(1)) stop "pull1"
    end if
    
    if(N_buf(2) > 0) then
      rc = f77_zmq_recv( zmq_socket_pull, double_buf, 8*N_buf(2), 0)
      if(rc /= 8*N_buf(2)) stop "pull2"
    end if
    
    if(N_buf(3) > 0) then
      rc = f77_zmq_recv( zmq_socket_pull, det_buf, 2*N_int*bit_kind*N_buf(3), 0)
      if(rc /= 2*N_int*bit_kind*N_buf(3)) stop "pull3"
    end if

    rc = f77_zmq_recv( zmq_socket_pull, task_id, 4, 0)
    if(rc /= 4) stop "pull4"
  end if
! Activate is zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  rc = f77_zmq_send( zmq_socket_pull, 'ok', 2, 0)
IRP_ENDIF

end subroutine
 
 
            
