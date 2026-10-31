! Fortran 90 counterpart of single_cell.cpp: S2M -> M2M -> M2L -> L2L -> L2P
! against direct summation. Same generator, so the same points and errors.
! Conventions: S2M(centre - source) M2M(parent - child) M2L(target - source)
!              L2L(child - parent) L2P(body - centre) P2P(target - source)
program single_cell
  use yukawa_f90
  implicit none
  integer, parameter :: ns = 20, nt = 8
  integer(8) :: state = 12345_8
  real(wp) :: cs(3), ct(3), off(3), skid(3,2), tkid(3,2)
  real(wp) :: src(3,ns), tgt(3,nt), q(ns), ref(4,nt), F(4), d(3), err, scale
  real(wp), allocatable :: M(:,:), Mp(:), L(:), Lk(:,:)
  integer :: sc(ns), tc(nt), i, j, k, p, n

  cs = [0.0_wp, 0.0_wp, 0.0_wp];  ct = [2.5_wp, 0.5_wp, -0.5_wp];  off = 0.1_wp
  skid(:,1) = cs + off;  skid(:,2) = cs - off
  tkid(:,1) = ct + off;  tkid(:,2) = ct - off

  do i = 1, ns
    sc(i) = mod(i - 1, 2) + 1
    do k = 1, 3
      src(k,i) = skid(k,sc(i)) + 0.1_wp * uniform()
    end do
    q(i) = uniform()
  end do
  do i = 1, nt
    tc(i) = mod(i - 1, 2) + 1
    do k = 1, 3
      tgt(k,i) = tkid(k,tc(i)) + 0.1_wp * uniform()
    end do
  end do

  ref = 0.0_wp
  do i = 1, nt
    do j = 1, ns
      d = tgt(:,i) - src(:,j)
      call P2P(d(1), d(2), d(3), q(j), ref(:,i))
    end do
  end do

  print '(a,i0,a,i0,a)', 'kernel exp(-r)/r, ', ns, ' sources, ', nt, ' targets'
  print '(a)', ' order  max relative error in potential'
  do p = FMMGEN_MINORDER, FMMGEN_MAXORDER, 2
    n = (p + 1) * (p + 2) * (p + 3) / 6
    allocate(M(n,2), Mp(n), L(n), Lk(n,2))
    M = 0.0_wp;  Mp = 0.0_wp;  L = 0.0_wp;  Lk = 0.0_wp
    do j = 1, ns
      d = skid(:,sc(j)) - src(:,j)
      call S2M(d(1), d(2), d(3), q(j), M(:,sc(j)), p)
    end do
    do k = 1, 2
      d = cs - skid(:,k)
      call M2M(d(1), d(2), d(3), M(:,k), Mp, p)
    end do
    d = ct - cs
    call M2L(d(1), d(2), d(3), Mp, L, p)
    do k = 1, 2
      d = tkid(:,k) - ct
      call L2L(d(1), d(2), d(3), L, Lk(:,k), p)
    end do
    err = 0.0_wp;  scale = 0.0_wp
    do i = 1, nt
      F = 0.0_wp
      d = tgt(:,i) - tkid(:,tc(i))
      call L2P(d(1), d(2), d(3), Lk(:,tc(i)), F, p)
      err = max(err, abs(F(1) - ref(1,i)))
      scale = max(scale, abs(ref(1,i)))
    end do
    print '(i6,2x,es10.3)', p, err / scale
    deallocate(M, Mp, L, Lk)
  end do

contains

  real(wp) function uniform()
    state = mod(state * 1103515245_8 + 12345_8, 2147483648_8)
    uniform = 2.0_wp * real(state, wp) / 2147483648.0_wp - 1.0_wp
  end function uniform

end program single_cell
