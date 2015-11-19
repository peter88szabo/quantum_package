subroutine get_excitation_degree(key1,key2,degree,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the excitation degree between two determinants
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key1(Nint,2)
  integer(bit_kind), intent(in)  :: key2(Nint,2)
  integer, intent(out)           :: degree
  
  integer                        :: l
  
  ASSERT (Nint > 0)
  
  degree = popcnt(xor( key1(1,1), key2(1,1))) +                      &
      popcnt(xor( key1(1,2), key2(1,2)))
  !DEC$ NOUNROLL
  do l=2,Nint
    degree = degree+ popcnt(xor( key1(l,1), key2(l,1))) +            &
        popcnt(xor( key1(l,2), key2(l,2)))
  enddo
  ASSERT (degree >= 0)
  degree = ishft(degree,-1)
  
end



subroutine get_excitation(det1,det2,exc,degree,phase,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the excitation operators between two determinants and the phase
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: det1(Nint,2)
  integer(bit_kind), intent(in)  :: det2(Nint,2)
  integer, intent(out)           :: exc(0:2,2,2)
  integer, intent(out)           :: degree
  double precision, intent(out)  :: phase
  ! exc(number,hole/particle,spin)
  ! ex :
  ! exc(0,1,1) = number of holes alpha
  ! exc(0,2,1) = number of particle alpha
  ! exc(0,2,2) = number of particle beta
  ! exc(1,2,1) = first particle alpha
  ! exc(1,1,1) = first hole     alpha
  ! exc(1,2,2) = first particle beta
  ! exc(1,1,2) = first hole     beta
  
  ASSERT (Nint > 0)
  
  !DIR$ FORCEINLINE
  call get_excitation_degree(det1,det2,degree,Nint)
  select case (degree)
      
    case (3:)
      degree = -1
      return
      
    case (2)
      call get_double_excitation(det1,det2,exc,phase,Nint)
      return
      
    case (1)
      call get_mono_excitation(det1,det2,exc,phase,Nint)
      return
      
    case(0)
      return
      
  end select
end

subroutine decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Decodes the exc arrays returned by get_excitation.
  ! h1,h2 : Holes
  ! p1,p2 : Particles
  ! s1,s2 : Spins (1:alpha, 2:beta)
  ! degree : Degree of excitation
  END_DOC
  integer, intent(in)            :: exc(0:2,2,2),degree
  integer, intent(out)           :: h1,h2,p1,p2,s1,s2
  ASSERT (degree > 0)
  ASSERT (degree < 3)
  
  select case(degree)
    case(2)
      if (exc(0,1,1) == 2) then
        h1 = exc(1,1,1)
        h2 = exc(2,1,1)
        p1 = exc(1,2,1)
        p2 = exc(2,2,1)
        s1 = 1
        s2 = 1
      else if (exc(0,1,2) == 2) then
        h1 = exc(1,1,2)
        h2 = exc(2,1,2)
        p1 = exc(1,2,2)
        p2 = exc(2,2,2)
        s1 = 2
        s2 = 2
      else
        h1 = exc(1,1,1)
        h2 = exc(1,1,2)
        p1 = exc(1,2,1)
        p2 = exc(1,2,2)
        s1 = 1
        s2 = 2
      endif
    case(1)
      if (exc(0,1,1) == 1) then
        h1 = exc(1,1,1)
        h2 = 0
        p1 = exc(1,2,1)
        p2 = 0
        s1 = 1
        s2 = 0
      else
        h1 = exc(1,1,2)
        h2 = 0
        p1 = exc(1,2,2)
        p2 = 0
        s1 = 2
        s2 = 0
      endif
    case(0)
      h1 = 0
      p1 = 0
      h2 = 0
      p2 = 0
      s1 = 0
      s2 = 0
  end select
end


subroutine get_double_excitation(det1,det2,exc,phase,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the two excitation operators between two doubly excited determinants and the phase
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: det1(Nint,2)
  integer(bit_kind), intent(in)  :: det2(Nint,2)
  integer, intent(out)           :: exc(0:2,2,2)
  double precision, intent(out)  :: phase
  integer                        :: tz
  integer                        :: l, ispin, idx_hole, idx_particle, ishift
  integer                        :: nperm
  integer                        :: i,j,k,m,n
  integer                        :: high, low
  integer                        :: a,b,c,d
  integer(bit_kind)              :: hole, particle, tmp
  double precision, parameter    :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)
  
  ASSERT (Nint > 0)
  nperm = 0
  exc(0,1,1) = 0
  exc(0,2,1) = 0
  exc(0,1,2) = 0
  exc(0,2,2) = 0
  do ispin = 1,2
    idx_particle = 0
    idx_hole = 0
    ishift = 1-bit_kind_size
    do l=1,Nint
      ishift = ishift + bit_kind_size
      if (det1(l,ispin) == det2(l,ispin)) then
        cycle
      endif
      tmp = xor( det1(l,ispin), det2(l,ispin) )
      particle = iand(tmp, det2(l,ispin))
      hole     = iand(tmp, det1(l,ispin))
      do while (particle /= 0_bit_kind)
        tz = trailz(particle)
        idx_particle = idx_particle + 1
        exc(0,2,ispin) = exc(0,2,ispin) + 1
        exc(idx_particle,2,ispin) = tz+ishift
        particle = iand(particle,particle-1_bit_kind)
      enddo
      if (iand(exc(0,1,ispin),exc(0,2,ispin))==2) then  ! exc(0,1,ispin)==2 or exc(0,2,ispin)==2
        exit
      endif
      do while (hole /= 0_bit_kind)
        tz = trailz(hole)
        idx_hole = idx_hole + 1
        exc(0,1,ispin) = exc(0,1,ispin) + 1
        exc(idx_hole,1,ispin) = tz+ishift
        hole = iand(hole,hole-1_bit_kind)
      enddo
      if (iand(exc(0,1,ispin),exc(0,2,ispin))==2) then ! exc(0,1,ispin)==2 or exc(0,2,ispin)
        exit
      endif
    enddo
    
    ! TODO : Voir si il faut sortir i,n,k,m du case.
    
    select case (exc(0,1,ispin))
      case(0)
        cycle
        
      case(1)
        low  = min(exc(1,1,ispin), exc(1,2,ispin))
        high = max(exc(1,1,ispin), exc(1,2,ispin))
        
        ASSERT (low > 0)
        j = ishft(low-1,-bit_kind_shift)+1   ! Find integer in array(Nint)
        n = iand(low,bit_kind_size-1)        ! mod(low,bit_kind_size)
        ASSERT (high > 0)
        k = ishft(high-1,-bit_kind_shift)+1
        m = iand(high,bit_kind_size-1)
        
        if (j==k) then
          nperm = nperm + popcnt(iand(det1(j,ispin),                 &
              iand( ibset(0_bit_kind,m-1)-1_bit_kind,                &
              ibclr(-1_bit_kind,n)+1_bit_kind ) ))
        else
          nperm = nperm + popcnt(iand(det1(k,ispin),                 &
              ibset(0_bit_kind,m-1)-1_bit_kind)) +                   &
              popcnt(iand(det1(j,ispin), ibclr(-1_bit_kind,n) +1_bit_kind))
          do i=j+1,k-1
            nperm = nperm + popcnt(det1(i,ispin))
          end do
        endif
        
      case (2)
        
        do i=1,2
          low  = min(exc(i,1,ispin), exc(i,2,ispin))
          high = max(exc(i,1,ispin), exc(i,2,ispin))
          
          ASSERT (low > 0)
          j = ishft(low-1,-bit_kind_shift)+1   ! Find integer in array(Nint)
          n = iand(low,bit_kind_size-1)        ! mod(low,bit_kind_size)
          ASSERT (high > 0)
          k = ishft(high-1,-bit_kind_shift)+1
          m = iand(high,bit_kind_size-1)
          
          if (j==k) then
            nperm = nperm + popcnt(iand(det1(j,ispin),               &
                iand( ibset(0_bit_kind,m-1)-1_bit_kind,              &
                ibclr(-1_bit_kind,n)+1_bit_kind ) ))
          else
            nperm = nperm + popcnt(iand(det1(k,ispin),               &
                ibset(0_bit_kind,m-1)-1_bit_kind)) +                 &
                popcnt(iand(det1(j,ispin), ibclr(-1_bit_kind,n) +1_bit_kind))
            do l=j+1,k-1
              nperm = nperm + popcnt(det1(l,ispin))
            end do
          endif
          
        enddo
        
        a = min(exc(1,1,ispin), exc(1,2,ispin))
        b = max(exc(1,1,ispin), exc(1,2,ispin))
        c = min(exc(2,1,ispin), exc(2,2,ispin))
        d = max(exc(2,1,ispin), exc(2,2,ispin))
        if (c>a .and. c<b .and. d>b) then
          nperm = nperm + 1
        endif
        exit
    end select
    
  enddo
  phase = phase_dble(iand(nperm,1))
  
end

subroutine get_mono_excitation(det1,det2,exc,phase,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the excitation operator between two singly excited determinants and the phase
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: det1(Nint,2)
  integer(bit_kind), intent(in)  :: det2(Nint,2)
  integer, intent(out)           :: exc(0:2,2,2)
  double precision, intent(out)  :: phase
  integer                        :: tz
  integer                        :: l, ispin, idx_hole, idx_particle, ishift
  integer                        :: nperm
  integer                        :: i,j,k,m,n
  integer                        :: high, low
  integer                        :: a,b,c,d
  integer(bit_kind)              :: hole, particle, tmp
  double precision, parameter    :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)
  
  ASSERT (Nint > 0)
  nperm = 0
  exc(0,1,1) = 0
  exc(0,2,1) = 0
  exc(0,1,2) = 0
  exc(0,2,2) = 0
  do ispin = 1,2
    ishift = 1-bit_kind_size
    do l=1,Nint
      ishift = ishift + bit_kind_size
      if (det1(l,ispin) == det2(l,ispin)) then
        cycle
      endif
      tmp = xor( det1(l,ispin), det2(l,ispin) )
      particle = iand(tmp, det2(l,ispin))
      hole     = iand(tmp, det1(l,ispin))
      if (particle /= 0_bit_kind) then
        tz = trailz(particle)
        exc(0,2,ispin) = 1
        exc(1,2,ispin) = tz+ishift
      endif
      if (hole /= 0_bit_kind) then
        tz = trailz(hole)
        exc(0,1,ispin) = 1
        exc(1,1,ispin) = tz+ishift
      endif
      
      if ( iand(exc(0,1,ispin),exc(0,2,ispin)) /= 1) then  ! exc(0,1,ispin)/=1 and exc(0,2,ispin) /= 1
        cycle
      endif
      
      low = min(exc(1,1,ispin),exc(1,2,ispin))
      high = max(exc(1,1,ispin),exc(1,2,ispin))
      
      ASSERT (low > 0)
      j = ishft(low-1,-bit_kind_shift)+1   ! Find integer in array(Nint)
      n = iand(low,bit_kind_size-1)        ! mod(low,bit_kind_size)
      ASSERT (high > 0)
      k = ishft(high-1,-bit_kind_shift)+1
      m = iand(high,bit_kind_size-1)
      if (j==k) then
        nperm = popcnt(iand(det1(j,ispin),                           &
            iand(ibset(0_bit_kind,m-1)-1_bit_kind,ibclr(-1_bit_kind,n)+1_bit_kind)))
      else
        nperm = nperm + popcnt(iand(det1(k,ispin),ibset(0_bit_kind,m-1)-1_bit_kind)) +&
            popcnt(iand(det1(j,ispin),ibclr(-1_bit_kind,n)+1_bit_kind))
        do i=j+1,k-1
          nperm = nperm + popcnt(det1(i,ispin))
        end do
      endif
      phase = phase_dble(iand(nperm,1))
      return
      
    enddo
  enddo
end





subroutine i_H_j(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij
  
  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: get_mo_bielec_integral
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem, phase,phase_2
  integer                        :: n_occ_alpha, n_occ_beta
  logical                        :: has_mipi(Nint*bit_kind_size)
  double precision               :: mipi(Nint*bit_kind_size), miip(Nint*bit_kind_size)
  PROVIDE mo_bielec_integrals_in_map mo_integrals_map
  
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)
  
  hij = 0.d0
  !DEC$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  select case (degree)
    case (2)
      call get_double_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha, mono beta
        hij = phase*get_mo_bielec_integral(                          &
            exc(1,1,1),                                              &
            exc(1,1,2),                                              &
            exc(1,2,1),                                              &
            exc(1,2,2) ,mo_integrals_map)
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij = phase*(get_mo_bielec_integral(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_map) -                          &
            get_mo_bielec_integral(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_map) )
      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij = phase*(get_mo_bielec_integral(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_map) -                          &
            get_mo_bielec_integral(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_map) )
      endif
    case (1)
      call get_mono_excitation(key_i,key_j,exc,phase,Nint)
      call bitstring_to_list(key_i(1,1), occ(1,1), n_occ_alpha, Nint)
      call bitstring_to_list(key_i(1,2), occ(1,2), n_occ_beta, Nint)
      has_mipi = .False.
      if (exc(0,1,1) == 1) then
        ! Mono alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        do k = 1, elec_alpha_num
          i = occ(k,1)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            miip(i) = get_mo_bielec_integral(m,i,i,p,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        do k = 1, elec_beta_num
          i = occ(k,2)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        
        do k = 1, elec_alpha_num
          hij = hij + mipi(occ(k,1)) - miip(occ(k,1))
        enddo
        do k = 1, elec_beta_num
          hij = hij + mipi(occ(k,2))
        enddo
        
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        do k = 1, elec_beta_num
          i = occ(k,2)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            miip(i) = get_mo_bielec_integral(m,i,i,p,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        do k = 1, elec_alpha_num
          i = occ(k,1)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        
        do k = 1, elec_alpha_num
          hij = hij + mipi(occ(k,1))
        enddo
        do k = 1, elec_beta_num
          hij = hij + mipi(occ(k,2)) - miip(occ(k,2))
        enddo
        
      endif
      hij = phase*(hij + mo_mono_elec_integral(m,p))
      
    case (0)
      hij = diag_H_mat_elem(key_i,Nint)
  end select
end



subroutine i_H_j_phase_out(key_i,key_j,Nint,hij,phase,exc,degree)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij, phase

  integer,intent(out)            :: exc(0:2,2,2)
  integer,intent(out)            :: degree
  double precision               :: get_mo_bielec_integral
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem
  integer                        :: n_occ_alpha, n_occ_beta
  logical                        :: has_mipi(Nint*bit_kind_size)
  double precision               :: mipi(Nint*bit_kind_size), miip(Nint*bit_kind_size)
  PROVIDE mo_bielec_integrals_in_map mo_integrals_map

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)

  hij = 0.d0
  !DEC$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  select case (degree)
    case (2)
      call get_double_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha, mono beta
        hij = phase*get_mo_bielec_integral(                          &
            exc(1,1,1),                                              &
            exc(1,1,2),                                              &
            exc(1,2,1),                                              &
            exc(1,2,2) ,mo_integrals_map)
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij = phase*(get_mo_bielec_integral(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_map) -                          &
            get_mo_bielec_integral(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_map) )
      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij = phase*(get_mo_bielec_integral(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_map) -                          &
            get_mo_bielec_integral(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_map) )
      endif
    case (1)
      call get_mono_excitation(key_i,key_j,exc,phase,Nint)
      call bitstring_to_list(key_i(1,1), occ(1,1), n_occ_alpha, Nint)
      call bitstring_to_list(key_i(1,2), occ(1,2), n_occ_beta, Nint)
      has_mipi = .False.
      if (exc(0,1,1) == 1) then
        ! Mono alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        do k = 1, elec_alpha_num
          i = occ(k,1)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            miip(i) = get_mo_bielec_integral(m,i,i,p,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        do k = 1, elec_beta_num
          i = occ(k,2)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo

        do k = 1, elec_alpha_num
          hij = hij + mipi(occ(k,1)) - miip(occ(k,1))
        enddo
        do k = 1, elec_beta_num
          hij = hij + mipi(occ(k,2))
        enddo

      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        do k = 1, elec_beta_num
          i = occ(k,2)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            miip(i) = get_mo_bielec_integral(m,i,i,p,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        do k = 1, elec_alpha_num
          i = occ(k,1)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo

        do k = 1, elec_alpha_num
          hij = hij + mipi(occ(k,1))
        enddo
        do k = 1, elec_beta_num
          hij = hij + mipi(occ(k,2)) - miip(occ(k,2))
        enddo

      endif
      hij = phase*(hij + mo_mono_elec_integral(m,p))

    case (0)
      hij = diag_H_mat_elem(key_i,Nint)
  end select
end



subroutine i_H_j_verbose(key_i,key_j,Nint,hij,hmono,hdouble)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij,hmono,hdouble
  
  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: get_mo_bielec_integral
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem, phase,phase_2
  integer                        :: n_occ_alpha, n_occ_beta
  logical                        :: has_mipi(Nint*bit_kind_size)
  double precision               :: mipi(Nint*bit_kind_size), miip(Nint*bit_kind_size)
  PROVIDE mo_bielec_integrals_in_map mo_integrals_map
  
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)
  
  hij = 0.d0
  hmono = 0.d0
  hdouble = 0.d0
  !DEC$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  select case (degree)
    case (2)
      call get_double_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha, mono beta
        hij = phase*get_mo_bielec_integral(                          &
            exc(1,1,1),                                              &
            exc(1,1,2),                                              &
            exc(1,2,1),                                              &
            exc(1,2,2) ,mo_integrals_map)
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij = phase*(get_mo_bielec_integral(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_map) -                          &
            get_mo_bielec_integral(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_map) )
      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij = phase*(get_mo_bielec_integral(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_map) -                          &
            get_mo_bielec_integral(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_map) )
      endif
    case (1)
      call get_mono_excitation(key_i,key_j,exc,phase,Nint)
      call bitstring_to_list(key_i(1,1), occ(1,1), n_occ_alpha, Nint)
      call bitstring_to_list(key_i(1,2), occ(1,2), n_occ_beta, Nint)
      has_mipi = .False.
      if (exc(0,1,1) == 1) then
        ! Mono alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        do k = 1, elec_alpha_num
          i = occ(k,1)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            miip(i) = get_mo_bielec_integral(m,i,i,p,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        do k = 1, elec_beta_num
          i = occ(k,2)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        
        do k = 1, elec_alpha_num
          hdouble = hdouble + mipi(occ(k,1)) - miip(occ(k,1))
        enddo
        do k = 1, elec_beta_num
          hdouble = hdouble + mipi(occ(k,2))
        enddo
        
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        do k = 1, elec_beta_num
          i = occ(k,2)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            miip(i) = get_mo_bielec_integral(m,i,i,p,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        do k = 1, elec_alpha_num
          i = occ(k,1)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral(m,i,p,i,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        
        do k = 1, elec_alpha_num
          hdouble = hdouble + mipi(occ(k,1))
        enddo
        do k = 1, elec_beta_num
          hdouble = hdouble + mipi(occ(k,2)) - miip(occ(k,2))
        enddo
        
      endif
      hmono = mo_mono_elec_integral(m,p)
      hij = phase*(hdouble + hmono)
      
    case (0)
      hij = diag_H_mat_elem(key_i,Nint)
  end select
end


subroutine create_minilist(key_mask, fullList, miniList, idx_miniList, N_fullList, N_miniList, Nint)
  use bitmasks
  implicit none
  
  integer(bit_kind), intent(in)            :: fullList(Nint, 2, N_fullList)
  integer, intent(in)                      :: N_fullList
  integer(bit_kind),intent(out)            :: miniList(Nint, 2, N_fullList)
  integer,intent(out)                      :: idx_miniList(N_fullList), N_miniList
  integer, intent(in)                      :: Nint
  integer(bit_kind)                        :: key_mask(Nint, 2)
  integer                                  :: ni, i, n_a, n_b, e_a, e_b
  
  
  n_a = 0
  n_b = 0
  do ni=1,nint
    n_a = n_a + popcnt(key_mask(ni,1))
    n_b = n_b + popcnt(key_mask(ni,2))
  end do
  
  if(n_a == 0) then
    N_miniList = N_fullList
    miniList(:,:,:) = fullList(:,:,:)
    do i=1,N_fullList
      idx_miniList(i) = i
    end do
    return
  end if
  
  N_miniList = 0
  
  do i=1,N_fullList
    e_a = n_a
    e_b = n_b
    do ni=1,nint
      e_a -= popcnt(iand(fullList(ni, 1, i), key_mask(ni, 1)))
      e_b -= popcnt(iand(fullList(ni, 2, i), key_mask(ni, 2)))
    end do
    
    if(e_a + e_b <= 2) then
      N_miniList = N_miniList + 1
      miniList(:,:,N_miniList) = fullList(:,:,i)
      idx_miniList(N_miniList) = i
    end if
  end do
end subroutine

subroutine i_H_psi_nominilist(key,keys,coef,Nint,Ndet,Ndet_max,Nstate,i_H_psi_array)
  use bitmasks
  implicit none
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate)
  double precision, intent(out)  :: i_H_psi_array(Nstate)
  
  integer                        :: i, ii,j
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: hij
  integer                        :: idx(0:Ndet)
  BEGIN_DOC
  ! <key|H|psi> for the various Nstates
  END_DOC
  
  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  i_H_psi_array = 0.d0
  
  call filter_connected_i_H_psi0(keys,key,Nint,Ndet,idx)
  do ii=1,idx(0)
    i = idx(ii)
    !DEC$ FORCEINLINE
    call i_H_j(keys(1,1,i),key,Nint,hij)
    do j = 1, Nstate
      i_H_psi_array(j) = i_H_psi_array(j) + coef(i,j)*hij
    enddo
  enddo
end


subroutine i_H_psi(key,keys,idx_key,N_minilist,coef,Nint,Ndet,Ndet_max,Nstate,i_H_psi_array)
  use bitmasks
  implicit none
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate,idx_key(Ndet), N_minilist
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate)
  double precision, intent(out)  :: i_H_psi_array(Nstate)
  
  integer                        :: i, ii,j, i_in_key, i_in_coef
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: hij
  integer                        :: idx(0:Ndet)
  BEGIN_DOC
  ! <key|H|psi> for the various Nstates
  END_DOC
  
  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  i_H_psi_array = 0.d0
  
  !call filter_connected_i_H_psi0(keys,key,Nint,Ndet,idx)
  call filter_connected_i_H_psi0(keys,key,Nint,N_minilist,idx)
  do ii=1,idx(0)
    !i = idx_key(idx(ii))
    i_in_key = idx(ii)
    i_in_coef = idx_key(idx(ii))
    !DEC$ FORCEINLINE
! !     call i_H_j(keys(1,1,i),key,Nint,hij)
! !     do j = 1, Nstate
! !       i_H_psi_array(j) = i_H_psi_array(j) + coef(i,j)*hij
! !     enddo
    call i_H_j(keys(1,1,i_in_key),key,Nint,hij)
    do j = 1, Nstate
      i_H_psi_array(j) = i_H_psi_array(j) + coef(i_in_coef,j)*hij
    enddo
  enddo
end

subroutine i_H_psi_sec_ord(key,keys,coef,Nint,Ndet,Ndet_max,Nstate,i_H_psi_array,idx_interaction,interactions)
  use bitmasks
  implicit none
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate)
  double precision, intent(out)  :: i_H_psi_array(Nstate)
  double precision, intent(out)  :: interactions(Ndet)
  integer,intent(out)            :: idx_interaction(0:Ndet)
  
  integer                        :: i, ii,j
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: hij
  integer                        :: idx(0:Ndet),n_interact
  BEGIN_DOC
  ! <key|H|psi> for the various Nstates
  END_DOC
  
  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  i_H_psi_array = 0.d0
  call filter_connected_i_H_psi0(keys,key,Nint,Ndet,idx)
  n_interact = 0
  do ii=1,idx(0)
    i = idx(ii)
    !DEC$ FORCEINLINE
    call i_H_j(keys(1,1,i),key,Nint,hij)
    if(dabs(hij).ge.1.d-8)then
     if(i.ne.1)then
      n_interact += 1
      interactions(n_interact) = hij
      idx_interaction(n_interact) = i
     endif
    endif
    do j = 1, Nstate
      i_H_psi_array(j) = i_H_psi_array(j) + coef(i,j)*hij
    enddo
  enddo
  idx_interaction(0) = n_interact
end


subroutine i_H_psi_SC2(key,keys,coef,Nint,Ndet,Ndet_max,Nstate,i_H_psi_array,idx_repeat)
  use bitmasks
  BEGIN_DOC
  ! <key|H|psi> for the various Nstate
  !
  ! returns in addition
  !
  ! the array of the index of the non connected determinants to key1
  !
  ! in order to know what double excitation can be repeated on key1
  !
  ! idx_repeat(0) is the number of determinants that can be used
  !
  ! to repeat the excitations
  END_DOC
  implicit none
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate)
  double precision, intent(out)  :: i_H_psi_array(Nstate)
  integer         , intent(out)  :: idx_repeat(0:Ndet)
  
  integer                        :: i, ii,j
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: hij
  integer                        :: idx(0:Ndet)
  
  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  i_H_psi_array = 0.d0
  call filter_connected_i_H_psi0_SC2(keys,key,Nint,Ndet,idx,idx_repeat)
  do ii=1,idx(0)
    i = idx(ii)
    !DEC$ FORCEINLINE
    call i_H_j(keys(1,1,i),key,Nint,hij)
    do j = 1, Nstate
      i_H_psi_array(j) = i_H_psi_array(j) + coef(i,j)*hij
    enddo
  enddo
end


subroutine i_H_psi_SC2_verbose(key,keys,coef,Nint,Ndet,Ndet_max,Nstate,i_H_psi_array,idx_repeat)
  use bitmasks
  BEGIN_DOC
  ! <key|H|psi> for the various Nstate
  !
  ! returns in addition
  !
  ! the array of the index of the non connected determinants to key1
  !
  ! in order to know what double excitation can be repeated on key1
  !
  ! idx_repeat(0) is the number of determinants that can be used
  !
  ! to repeat the excitations
  END_DOC
  implicit none
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate)
  double precision, intent(out)  :: i_H_psi_array(Nstate)
  integer         , intent(out)  :: idx_repeat(0:Ndet)
  
  integer                        :: i, ii,j
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: hij
  integer                        :: idx(0:Ndet)
  
  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  i_H_psi_array = 0.d0
  call filter_connected_i_H_psi0_SC2(keys,key,Nint,Ndet,idx,idx_repeat)
  print*,'--------'
  do ii=1,idx(0)
    print*,'--'
    i = idx(ii)
    !DEC$ FORCEINLINE
    call i_H_j(keys(1,1,i),key,Nint,hij)
    if (i==1)then
     print*,'i==1 !!'
    endif
    print*,coef(i,1) * hij,coef(i,1),hij
    do j = 1, Nstate
      i_H_psi_array(j) = i_H_psi_array(j) + coef(i,j)*hij
    enddo
    print*,i_H_psi_array(1)
  enddo
  print*,'------'
end



subroutine get_excitation_degree_vector(key1,key2,degree,Nint,sze,idx)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Applies get_excitation_degree to an array of determinants
  END_DOC
  integer, intent(in)            :: Nint, sze
  integer(bit_kind), intent(in)  :: key1(Nint,2,sze)
  integer(bit_kind), intent(in)  :: key2(Nint,2)
  integer, intent(out)           :: degree(sze)
  integer, intent(out)           :: idx(0:sze)
  
  integer                        :: i,l,d,m
  
  ASSERT (Nint > 0)
  ASSERT (sze > 0)
  
  l=1
  if (Nint==1) then
    
    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +       &
          popcnt(xor( key1(1,2,i), key2(1,2)))
      if (d > 4) then
        cycle
      else
        degree(l) = ishft(d,-1)
        idx(l) = i
        l = l+1
      endif
    enddo
    
  else if (Nint==2) then
    
    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +                     &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2)))
      if (d > 4) then
        cycle
      else
        degree(l) = ishft(d,-1)
        idx(l)    = i
        l         = l+1
      endif
    enddo
    
  else if (Nint==3) then
    
    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +                     &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2))) +                     &
          popcnt(xor( key1(3,1,i), key2(3,1))) +                     &
          popcnt(xor( key1(3,2,i), key2(3,2)))
      if (d > 4) then
        cycle
      else
        degree(l) = ishft(d,-1)
        idx(l)    = i
        l         = l+1
      endif
    enddo
    
  else
    
    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = 0
      !DEC$ LOOP COUNT MIN(4)
      do m=1,Nint
        d = d + popcnt(xor( key1(m,1,i), key2(m,1)))                 &
              + popcnt(xor( key1(m,2,i), key2(m,2)))
      enddo
      if (d > 4) then
        cycle
      else
        degree(l) = ishft(d,-1)
        idx(l)    = i
        l         = l+1
      endif
    enddo
    
  endif
  idx(0) = l-1
end




double precision function diag_H_mat_elem(det_in,Nint)
  implicit none
  BEGIN_DOC
  ! Computes <i|H|i>
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)
  
  integer(bit_kind)              :: hole(Nint,2)
  integer(bit_kind)              :: particle(Nint,2)
  integer                        :: i, nexc(2), ispin
  integer                        :: occ_particle(Nint*bit_kind_size,2)
  integer                        :: occ_hole(Nint*bit_kind_size,2)
  integer(bit_kind)              :: det_tmp(Nint,2)
  integer                        :: na, nb
  
  ASSERT (Nint > 0)
  ASSERT (sum(popcnt(det_in(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(det_in(:,2))) == elec_beta_num)
  
  nexc(1) = 0
  nexc(2) = 0
  do i=1,Nint
    hole(i,1)     =  xor(det_in(i,1),ref_bitmask(i,1))
    hole(i,2)     =  xor(det_in(i,2),ref_bitmask(i,2))
    particle(i,1) = iand(hole(i,1),det_in(i,1))
    particle(i,2) = iand(hole(i,2),det_in(i,2))
    hole(i,1)     = iand(hole(i,1),ref_bitmask(i,1))
    hole(i,2)     = iand(hole(i,2),ref_bitmask(i,2))
    nexc(1)      += popcnt(hole(i,1))
    nexc(2)      += popcnt(hole(i,2))
  enddo
  
  diag_H_mat_elem = ref_bitmask_energy
  if (nexc(1)+nexc(2) == 0) then
    return
  endif
  
  !call debug_det(det_in,Nint)
  integer                        :: tmp
  call bitstring_to_list(particle(1,1), occ_particle(1,1), tmp, Nint)
  ASSERT (tmp == nexc(1))
  call bitstring_to_list(particle(1,2), occ_particle(1,2), tmp, Nint)
  ASSERT (tmp == nexc(2))
  call bitstring_to_list(hole(1,1), occ_hole(1,1), tmp, Nint)
  ASSERT (tmp == nexc(1))
  call bitstring_to_list(hole(1,2), occ_hole(1,2), tmp, Nint)
  ASSERT (tmp == nexc(2))
  
  det_tmp = ref_bitmask
  do ispin=1,2
    na = elec_num_tab(ispin)
    nb = elec_num_tab(iand(ispin,1)+1)
    do i=1,nexc(ispin)
      !DIR$ FORCEINLINE
      call ac_operator( occ_particle(i,ispin), ispin, det_tmp, diag_H_mat_elem, Nint,na,nb)
      !DIR$ FORCEINLINE
      call a_operator ( occ_hole    (i,ispin), ispin, det_tmp, diag_H_mat_elem, Nint,na,nb)
    enddo
  enddo
end

subroutine a_operator(iorb,ispin,key,hjj,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Needed for diag_H_mat_elem
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hjj
  
  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: other_spin
  integer                        :: k,l,i
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k > 0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  key(k,ispin) = ibclr(key(k,ispin),l)
  other_spin = iand(ispin,1)+1
  
  !DIR$ FORCEINLINE
  call get_occ_from_key(key,occ,Nint)
  na -= 1
  
  hjj -= mo_mono_elec_integral(iorb,iorb)
  
  ! Same spin
  do i=1,na
    hjj -= mo_bielec_integral_jj_anti(occ(i,ispin),iorb)
  enddo
  
  ! Opposite spin
  do i=1,nb
    hjj -= mo_bielec_integral_jj(occ(i,other_spin),iorb)
  enddo
  
end


subroutine ac_operator(iorb,ispin,key,hjj,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Needed for diag_H_mat_elem
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hjj
  
  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: other_spin
  integer                        :: k,l,i
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  integer                        :: tmp
  !DIR$ FORCEINLINE
  call bitstring_to_list(key(1,1), occ(1,1), tmp, Nint)
  ASSERT (tmp == elec_alpha_num)
  !DIR$ FORCEINLINE
  call bitstring_to_list(key(1,2), occ(1,2), tmp, Nint)
  ASSERT (tmp == elec_beta_num)
  
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k > 0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  key(k,ispin) = ibset(key(k,ispin),l)
  other_spin = iand(ispin,1)+1
  
  hjj += mo_mono_elec_integral(iorb,iorb)
  
  ! Same spin
  do i=1,na
    hjj += mo_bielec_integral_jj_anti(occ(i,ispin),iorb)
  enddo
  
  ! Opposite spin
  do i=1,nb
    hjj += mo_bielec_integral_jj(occ(i,other_spin),iorb)
  enddo
  na += 1
end

subroutine get_occ_from_key(key,occ,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns a list of occupation numbers from a bitstring
  END_DOC
  integer(bit_kind), intent(in)  :: key(Nint,2)
  integer          , intent(in)  :: Nint
  integer         , intent(out)  :: occ(Nint*bit_kind_size,2)
  integer                        :: tmp
  
  call bitstring_to_list(key(1,1), occ(1,1), tmp, Nint)
  call bitstring_to_list(key(1,2), occ(1,2), tmp, Nint)
  
end

subroutine H_u_0(v_0,u_0,H_jj,n,keys_tmp,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0 = H|u_0>
  !
  ! n : number of determinants
  !
  ! H_jj : array of <j|H|j>
  END_DOC
  integer, intent(in)            :: n,Nint
  double precision, intent(out)  :: v_0(n)
  double precision, intent(in)   :: u_0(n)
  double precision, intent(in)   :: H_jj(n)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)
  integer, allocatable           :: idx(:)
  double precision               :: hij
  double precision, allocatable  :: vt(:)
  integer                        :: i,j,k,l, jj,ii
  integer                        :: i0, j0
  
  integer                        :: shortcut(0:n+1), sort_idx(n)
  integer(bit_kind)              :: sorted(Nint,n), version(Nint,n)
  
  
  integer                        :: sh, sh2, ni, exa, ext, org_i, org_j, endi
!   
  

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (n>0)
  PROVIDE ref_bitmask_energy
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP PRIVATE(i,hij,j,k,idx,jj,vt,ii,sh,sh2,ni,exa,ext,org_i,org_j,endi)                             &
      !$OMP SHARED(n,H_jj,u_0,keys_tmp,Nint,v_0,davidson_threshold,sorted,shortcut,sort_idx,version,davidson_criterion_is_built)
  allocate(idx(0:n), vt(n))
  Vt = 0.d0
  v_0 = 0.d0
  
  !$OMP SINGLE
  call sort_dets_ab_v(keys_tmp, sorted, sort_idx, shortcut, version, n, Nint)
  !$OMP END SINGLE
  
  !$OMP DO SCHEDULE(dynamic)
  do sh=1,shortcut(0)
  do sh2=1,sh
    exa = 0
    do ni=1,Nint
      exa += popcnt(xor(version(ni,sh), version(ni,sh2)))
    end do
    if(exa > 2) then
      cycle
    end if
    
    do i=shortcut(sh),shortcut(sh+1)-1
      if(sh==sh2) then
        endi = i-1
      else
        endi = shortcut(sh2+1)-1
      end if
      
      do j=shortcut(sh2),endi
       ext = exa
        do ni=1,Nint
          ext += popcnt(xor(sorted(ni,i), sorted(ni,j)))
        end do
        if(ext <= 4) then
          org_i = sort_idx(i)
          org_j = sort_idx(j)
          if ( dabs(u_0(org_j)) + dabs(u_0(org_i)) > davidson_threshold ) then
            call i_H_j(keys_tmp(1,1,org_j),keys_tmp(1,1,org_i),Nint,hij)
            vt (org_i) = vt (org_i) + hij*u_0(org_j)
            vt (org_j) = vt (org_j) + hij*u_0(org_i)
          endif
        endif
      enddo
    enddo
  enddo
  enddo
  !$OMP END DO
  
  !$OMP SINGLE
  call sort_dets_ba_v(keys_tmp, sorted, sort_idx, shortcut, version, n, Nint)
  !$OMP END SINGLE
 
  !$OMP DO SCHEDULE(dynamic)
  do sh=1,shortcut(0)
    do i=shortcut(sh),shortcut(sh+1)-1
      do j=shortcut(sh),i-1
        ext = 0
        do ni=1,Nint
          ext += popcnt(xor(sorted(ni,i), sorted(ni,j)))
        end do
        if(ext == 4) then
          org_i = sort_idx(i)
          org_j = sort_idx(j)
          if ( dabs(u_0(org_j)) + dabs(u_0(org_i)) > davidson_threshold ) then
            call i_H_j(keys_tmp(1,1,org_j),keys_tmp(1,1,org_i),Nint,hij)
            vt (org_i) = vt (org_i) + hij*u_0(org_j)
            vt (org_j) = vt (org_j) + hij*u_0(org_i)
          end if
        end if
      end do
    end do
  enddo
  !$OMP END DO
  
  !$OMP CRITICAL
  do i=1,n
    v_0(i) = v_0(i) + vt(i)
  enddo
  !$OMP END CRITICAL
  deallocate(idx,vt)
  !$OMP END PARALLEL
  do i=1,n
    v_0(i) += H_jj(i) * u_0(i)
  enddo
end
