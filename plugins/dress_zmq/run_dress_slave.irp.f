use bitmasks

BEGIN_PROVIDER [ integer, fragment_count ]
  implicit none
  BEGIN_DOC
    ! Number of fragments for the deterministic part
  END_DOC
  fragment_count = 1
END_PROVIDER


subroutine run_dress_slave(thread,iproc,energy)
  use f77_zmq
  implicit none

  double precision, intent(in)    :: energy(N_states_diag)
  integer,  intent(in)            :: thread, iproc
  integer                        :: rc, i, subset, i_generator(60)

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
  
  task(:) = CHAR(0)

  allocate(int_buf(N_dress_int_buffer)) 
  allocate(double_buf(N_dress_double_buffer)) 
  allocate(det_buf(N_int, 2, N_dress_det_buffer)) 
  allocate(delta_ij_loc(N_states,N_det,2)) 
  
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
  do
    call get_task_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id, task)
    if(task_id /= 0) then
      task = trim(task)//' 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
      
      i_generator = 0
      read (task,*) subset, i_generator
      if(i_generator(size(i_generator)) /= 0) stop "i_generator buffer too small"
      delta_ij_loc = 0d0
      i=1
      do while(i_generator(i) /= 0)
        call generator_start(i_generator(i), iproc)
        call alpha_callback(delta_ij_loc, i_generator(i), subset, iproc)
        call generator_done(i_generator(i), int_buf, double_buf, det_buf, N_buf, iproc)
        last = (i_generator(i+1) == 0)
        call push_dress_results(zmq_socket_push, i_generator(i), last, delta_ij_loc, int_buf, double_buf, det_buf, N_buf, task_id)
        i += 1
      end do
      call task_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id)
    else
      exit
    end if
  end do
  call disconnect_from_taskserver(zmq_to_qp_run_socket,zmq_socket_push,worker_id)
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_push_socket(zmq_socket_push,thread)
end subroutine


! BEGIN_PROVIDER [ integer, dress_int_buffer, (N_dress_int_buffer) ]
!&BEGIN_PROVIDER [ double precision, dress_double_buffer, (N_dress_double_buffer) ]
!&BEGIN_PROVIDER [ integer(bit_kind), dress_det_buffer, (N_int, 2, N_dress_det_buffer) ]
!  implicit none
!  
!  dress_int_buffer = 0
!  dress_double_buffer = 0d0
 ! dress_det_buffer = 0_bit_kind
!END_PROVIDER


!subroutine pull_dress_results(zmq_socket_pull, ind, delta_loc, int_buf, double_buf, det_buf, N_buf, task_id, felem)
subroutine push_dress_results(zmq_socket_push, ind, last, delta_loc, int_buf, double_buf, det_buf, N_bufi, task_id)
  use f77_zmq
  implicit none

  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  double precision, intent(in)   :: delta_loc(N_states, N_det, 2)
  double precision, intent(in) :: double_buf(*)
  logical, intent(in) :: last
  integer, intent(in) :: int_buf(*)
  integer(bit_kind), intent(in) :: det_buf(N_int, 2, *)
  integer, intent(in) :: N_bufi(3)
  integer :: N_buf(3)
  integer, intent(in) :: ind, task_id
  integer :: rc, i, j, felem
  double precision :: vast_emptiness(N_states)
  integer :: fillness

  vast_emptiness = 0d0
  felem = 1

  dloop : do i=1, N_det
    do j=1,N_states
      if(delta_loc(j,i,1) /= 0d0 .or. delta_loc(j,i,2) /= 0d0) then
        felem = i
        exit dloop
      end if
    end do
  end do dloop
  
  if(last) then
    fillness = 0
    do i=felem,N_det
      do j=1,N_states
        if(delta_loc(j,i,1) /= 0d0 .or. delta_loc(j,i,2) /= 0d0) then
          fillness += 1
        end if
      end do
    end do
    !print *, "FILLNESS", float(fillness) / float((N_det-felem+1)*N_states)
  end if


  rc = f77_zmq_send( zmq_socket_push, ind, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop "push"
  
  rc = f77_zmq_send( zmq_socket_push, last, 1, ZMQ_SNDMORE)
  if(rc /= 1) stop "push"
  
  if(last) then
    rc = f77_zmq_send( zmq_socket_push, felem, 4, ZMQ_SNDMORE)
    if(rc /= 4) stop "push"
 
    rc = f77_zmq_send( zmq_socket_push, delta_loc(1,felem,1), 8*N_states*(N_det-felem+1), ZMQ_SNDMORE)
    if(rc /= 8*N_states*(N_det+1-felem)) stop "push"

    rc = f77_zmq_send( zmq_socket_push, delta_loc(1,felem,2), 8*N_states*(N_det-felem+1), ZMQ_SNDMORE)
    if(rc /= 8*N_states*(N_det+1-felem)) stop "push"
  else
    rc = f77_zmq_send( zmq_socket_push, N_det, 4, ZMQ_SNDMORE)
    if(rc /= 4) stop "push"
 
    rc = f77_zmq_send( zmq_socket_push, vast_emptiness, 8*N_states, ZMQ_SNDMORE)
    if(rc /= 8*N_states) stop "push"

    rc = f77_zmq_send( zmq_socket_push, vast_emptiness, 8*N_states, ZMQ_SNDMORE)
    if(rc /= 8*N_states) stop "push" 
  end if
 
 
  N_buf = N_bufi

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

! Activate is zmq_socket_push is a REQ
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(2) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
IRP_ENDIF

end subroutine


subroutine pull_dress_results(zmq_socket_pull, ind, last, delta_loc, int_buf, double_buf, det_buf, N_buf, task_id, felem)
  use f77_zmq
  implicit none
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  logical, intent(out) :: last
  double precision, intent(inout) :: delta_loc(N_states, N_det, 2)
  double precision, intent(out) :: double_buf(*)
  integer, intent(out) :: int_buf(*)
  integer(bit_kind), intent(out) :: det_buf(N_int, 2, *)
  integer, intent(out) :: felem
  integer, intent(out) :: ind
  integer, intent(out) :: task_id
  integer :: rc, i
  integer, intent(out) :: N_buf(3)


  rc = f77_zmq_recv( zmq_socket_pull, ind, 4, 0)
  if(rc /= 4) stop "pulla"
 
  rc = f77_zmq_recv( zmq_socket_pull, last, 1, 0)
  if(rc /= 1) stop "pulla"
  
  rc = f77_zmq_recv( zmq_socket_pull, felem, 4, 0)
  if(rc /= 4) stop "pullb"
  
  delta_loc(:,:felem,:) = 0d0

  rc = f77_zmq_recv( zmq_socket_pull, delta_loc(1,felem,1), N_states*8*(N_det+1-felem), 0)
  if(rc /= 8*N_states*(N_det+1-felem)) stop "pullc"
  
  rc = f77_zmq_recv( zmq_socket_pull, delta_loc(1,felem,2), N_states*8*(N_det+1-felem), 0)
  if(rc /= 8*N_states*(N_det+1-felem)) stop "pulld"


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

! Activate is zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  rc = f77_zmq_send( zmq_socket_pull, 'ok', 2, 0)
IRP_ENDIF

end subroutine
 
 
            
