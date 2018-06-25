subroutine run_selection_slave(thread,iproc,energy)
  call run_selection_slave_new(thread,iproc,energy)
end

subroutine run_selection_slave_new(thread,iproc,energy)
  use f77_zmq
  use selection_types
  implicit none

  double precision, intent(in)    :: energy(N_states_diag)
  integer,  intent(in)            :: thread, iproc
  integer                        :: rc, i, N
  logical                        :: buffer_ready

  integer                        :: worker_id, ltask
  character*(512), allocatable   :: task(:)
  integer, allocatable           :: task_id(:)

  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_push_socket
  integer(ZMQ_PTR)               :: zmq_socket_push

  type(selection_buffer) :: buf, buf2
  logical :: done

  double precision,allocatable :: pt2(:,:)
  integer :: n_tasks, k, n_tasks_max
  integer, allocatable :: i_generator(:), subset(:)

  PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_bilinear_matrix_rows psi_det_sorted_order psi_bilinear_matrix_order
  PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  PROVIDE psi_bilinear_matrix_transp_order

  buffer_ready = .False.
  n_tasks_max = N_det_generators/100+1
  allocate(task_id(n_tasks_max), task(n_tasks_max))
  allocate(pt2(N_states,n_tasks_max), i_generator(n_tasks_max), subset(n_tasks_max))

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()

  integer, external :: connect_to_taskserver
  if (connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread) == -1) then
    call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
    return
  endif

  zmq_socket_push      = new_zmq_push_socket(thread)

  buf%N = 0
  n_tasks = 1
  call create_selection_buffer(0, 0, buf)
  done = .False.
  do while (.not.done)

    n_tasks = max(1,n_tasks)
    n_tasks = min(n_tasks,n_tasks_max)

    integer, external :: get_tasks_from_taskserver
    if (get_tasks_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id, task, n_tasks) == -1) then
      exit
    endif
    done = task_id(n_tasks) == 0
    if (done) n_tasks = n_tasks-1
    if (n_tasks == 0) exit

    do k=1,n_tasks
      read (task(k),*) subset(k), i_generator(k), N
    enddo

    if(buf%N == 0) then
        ! Only first time 
        call create_selection_buffer(N, N*2, buf)
        call create_selection_buffer(N, N*2, buf2)
        buffer_ready = .True.
    endif

    double precision :: time0, time1
    call wall_time(time0)
    do k=1,n_tasks
        pt2(:,k) = 0.d0
        buf%cur = 0
        call select_connected(i_generator(k),energy,pt2(1,k),buf,subset(k))
    enddo
    call wall_time(time1)

    integer, external :: tasks_done_to_taskserver
    if (tasks_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id,n_tasks) == -1) then
      done = .true.
    endif
    call sort_selection_buffer(buf)
    call merge_selection_buffers(buf,buf2)
    call push_selection_results(zmq_socket_push, pt2, buf, task_id, n_tasks)
    buf%mini = buf2%mini
    pt2(:,:) = 0d0
    buf%cur = 0

!    ! Try to adjust n_tasks around 5 second per job
!    n_tasks = min(n_tasks,int( 5.d0 * dble(n_tasks) / (time1 - time0 + 1.d-9)))+1
     n_tasks = n_tasks+1
  end do

  integer, external :: disconnect_from_taskserver
  if (disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id) == -1) then
    continue
  endif

  call end_zmq_push_socket(zmq_socket_push,thread)
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call delete_selection_buffer(buf)

end

subroutine run_selection_slave_old(thread,iproc,energy)
  use f77_zmq
  use selection_types
  implicit none

  double precision, intent(in)    :: energy(N_states)
  integer,  intent(in)            :: thread, iproc
  integer                        :: rc, i

  integer                        :: worker_id, task_id(1), ctask, ltask
  character*(512)                :: task

  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_push_socket
  integer(ZMQ_PTR)               :: zmq_socket_push

  type(selection_buffer) :: buf, buf2
  logical :: done, buffer_ready
  double precision :: pt2(N_states)

  PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_bilinear_matrix_rows psi_det_sorted_order psi_bilinear_matrix_order
  PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  PROVIDE psi_bilinear_matrix_transp_order

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()

  integer, external :: connect_to_taskserver 
  if (connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread) == -1) then 
    call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket) 
    return
  endif

  zmq_socket_push      = new_zmq_push_socket(thread)

  buf%N = 0
  buffer_ready = .False.
  ctask = 1
  pt2(:) = 0d0

  do
    integer, external :: get_task_from_taskserver
    if (get_task_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id(ctask), task) == -1) then
      exit
    endif
    done = task_id(ctask) == 0
    if (done) then
      ctask = ctask - 1
    else
      integer :: i_generator, N, subset
      read(task,*) subset, i_generator, N
      if(buf%N == 0) then
        ! Only first time 
        call create_selection_buffer(N, N*2, buf)
        call create_selection_buffer(N, N*2, buf2)
        buffer_ready = .True.
      else
        ASSERT (N == buf%N)
      end if
      call select_connected(i_generator,energy,pt2,buf,subset)
    endif

    integer, external :: task_done_to_taskserver

    if(done .or. ctask == size(task_id)) then
      do i=1, ctask
         if (task_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id(i)) == -1) then
           call sleep(1)
          if (task_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id(i)) == -1) then
            ctask = 0
            done = .true.
            exit
          endif
         endif
      end do
      if(ctask > 0) then
        call sort_selection_buffer(buf)
        call merge_selection_buffers(buf,buf2)
        call push_selection_results(zmq_socket_push, pt2, buf, task_id(1), ctask)
        buf%mini = buf2%mini
        pt2(:) = 0d0
        buf%cur = 0
      end if
      ctask = 0
    end if

    if(done) exit
    ctask = ctask + 1
  end do


  integer, external :: disconnect_from_taskserver 
  if (disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id) == -1) then 
    continue 
  endif 

  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_push_socket(zmq_socket_push,thread)
  if (buffer_ready) then
    call delete_selection_buffer(buf)
    call delete_selection_buffer(buf2)
  endif
end subroutine


subroutine push_selection_results(zmq_socket_push, pt2, b, task_id, ntask)
  use f77_zmq
  use selection_types
  implicit none

  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  double precision, intent(in)   :: pt2(N_states)
  type(selection_buffer), intent(inout) :: b
  integer, intent(in) :: ntask, task_id(*)
  integer :: rc

  rc = f77_zmq_send( zmq_socket_push, b%cur, 4, ZMQ_SNDMORE)
  if(rc /= 4) then
    print *,  'f77_zmq_send( zmq_socket_push, b%cur, 4, ZMQ_SNDMORE)'
  endif

  if (b%cur > 0) then

      rc = f77_zmq_send( zmq_socket_push, pt2, 8*N_states, ZMQ_SNDMORE)
      if(rc /= 8*N_states) then
        print *,  'f77_zmq_send( zmq_socket_push, pt2, 8*N_states, ZMQ_SNDMORE)'
      endif

      rc = f77_zmq_send( zmq_socket_push, b%val(1), 8*b%cur, ZMQ_SNDMORE)
      if(rc /= 8*b%cur) then
        print *,  'f77_zmq_send( zmq_socket_push, b%val(1), 8*b%cur, ZMQ_SNDMORE)'
      endif

      rc = f77_zmq_send( zmq_socket_push, b%det(1,1,1), bit_kind*N_int*2*b%cur, ZMQ_SNDMORE)
      if(rc /= bit_kind*N_int*2*b%cur) then
        print *,  'f77_zmq_send( zmq_socket_push, b%det(1,1,1), bit_kind*N_int*2*b%cur, ZMQ_SNDMORE)'
      endif

  endif

  rc = f77_zmq_send( zmq_socket_push, ntask, 4, ZMQ_SNDMORE)
  if(rc /= 4) then
    print *,  'f77_zmq_send( zmq_socket_push, ntask, 4, ZMQ_SNDMORE)'
  endif

  rc = f77_zmq_send( zmq_socket_push, task_id(1), ntask*4, 0)
  if(rc /= 4*ntask) then
    print *,  'f77_zmq_send( zmq_socket_push, task_id(1), ntask*4, 0)'
  endif

! Activate is zmq_socket_push is a REQ
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(2) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
  if ((rc /= 2).and.(ok(1:2) /= 'ok')) then
    print *,  irp_here//': error in receiving ok'
    stop -1
  endif
IRP_ENDIF

end subroutine


subroutine pull_selection_results(zmq_socket_pull, pt2, val, det, N, task_id, ntask)
  use f77_zmq
  use selection_types
  implicit none
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  double precision, intent(inout) :: pt2(N_states)
  double precision, intent(out) :: val(*)
  integer(bit_kind), intent(out) :: det(N_int, 2, *)
  integer, intent(out) :: N, ntask, task_id(*)
  integer :: rc, rn, i

  rc = f77_zmq_recv( zmq_socket_pull, N, 4, 0)
  if(rc /= 4) then
    print *,  'f77_zmq_recv( zmq_socket_pull, N, 4, 0)'
  endif

  if (N>0) then
      rc = f77_zmq_recv( zmq_socket_pull, pt2, N_states*8, 0)
      if(rc /= 8*N_states) then
        print *,  'f77_zmq_recv( zmq_socket_pull, pt2, N_states*8, 0)'
      endif

      rc = f77_zmq_recv( zmq_socket_pull, val(1), 8*N, 0)
      if(rc /= 8*N) then
        print *,  'f77_zmq_recv( zmq_socket_pull, val(1), 8*N, 0)'
      endif

      rc = f77_zmq_recv( zmq_socket_pull, det(1,1,1), bit_kind*N_int*2*N, 0)
      if(rc /= bit_kind*N_int*2*N) then
        print *,  'f77_zmq_recv( zmq_socket_pull, det(1,1,1), bit_kind*N_int*2*N, 0)'
      endif
  else
      pt2(:) = 0.d0
  endif

  rc = f77_zmq_recv( zmq_socket_pull, ntask, 4, 0)
  if(rc /= 4) then
    print *,  'f77_zmq_recv( zmq_socket_pull, ntask, 4, 0)'
  endif

  rc = f77_zmq_recv( zmq_socket_pull, task_id(1), ntask*4, 0)
  if(rc /= 4*ntask) then
    print *,  'f77_zmq_recv( zmq_socket_pull, task_id(1), ntask*4, 0)'
  endif

! Activate is zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  rc = f77_zmq_send( zmq_socket_pull, 'ok', 2, 0)
  if (rc /= 2) then
    print *,  irp_here//': error in sending ok'
    stop -1
  endif
IRP_ENDIF
end subroutine
 
 

