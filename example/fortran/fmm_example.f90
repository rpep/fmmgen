! Fortran 90 counterpart of the C++ example in ../: a Cartesian FMM and
! Barnes-Hut tree code built on the generated operators in operators.f90.
!
! Same algorithms, command line and output files as the C++ program:
!   - D = 3 (octree) with the full or the harmonic-compressed operators, or
!     D = 2 (quadtree, particles in z = 0) with the planar operators;
!   - FMM (dual-tree traversal, M2L + P2P) and Barnes-Hut (M2P + P2P);
!   - OpenMP over exclusive targets, so results do not depend on the number
!     of threads.
! Particle positions come from random_number, so errors differ slightly from
! the C++ program's.
!
! Conventions (copied from the C++ driver):
!   S2M(centre - source)   M2M(parent - child)   M2L(target - source)
!   L2L(child - parent)    L2P(body - centre)    M2P(body - cell)
!
! operators.f90 must have been generated with compress=.True. and
! planar=.True. (as example.py does), since the variant is chosen at run time.
module fmm_tree
!$ use omp_lib
  use operators
  implicit none
  private
  public :: build_tree, solve_fmm, solve_bh, solve_direct, refresh_sources
  public :: select_variant, variant_name, tree_ncells

  integer, parameter :: osz = FMMGEN_OUTPUTSIZE
  integer, parameter :: ssz = FMMGEN_SOURCESIZE

  integer, parameter :: V_FULL = 1, V_COMPRESSED = 2, V_PLANAR = 3
  integer :: variant = V_FULL

  integer :: dim = 3, nch = 8
  integer :: n = 0, ncrit = 0, order = 0, ncell = 0, ms = 0, ls = 0
  real(wp) :: theta = 0.0d0

  ! Cells. A cell is a leaf when nleaf < ncrit, as in the C++ code.
  real(wp), allocatable :: cen(:,:), rad(:), rmx(:)
  integer, allocatable :: nleaf(:), nchild(:), child(:,:), boff(:), parent(:), lvl(:)
  integer :: maxlev = 0
  ! Cells grouped by level: level l is lev_cells(lev_start(l) : lev_start(l+1)-1).
  integer, allocatable :: lev_start(:), lev_cells(:)

  ! Particles in tree (Morton) order. bz is all zero for D = 2.
  real(wp), allocatable :: bx(:), by(:), bz(:), bS(:), Fm(:)
  integer, allocatable :: perm(:)

  real(wp), allocatable :: Mc(:,:), Lc(:,:)

  ! Interaction lists, grouped by target cell: the sources of target A are
  ! src(start(A)) .. src(start(A+1)-1).
  integer, allocatable :: m2l_start(:), m2l_src(:), p2p_start(:), p2p_src(:)

  ! Scratch while building.
  integer, allocatable :: idx(:), tmp(:)
  real(wp), allocatable :: pos(:,:)
  integer, allocatable :: il_t(:), il_s(:), il_k(:)

  ! A growable list of interactions (kind 1 = M2L, 2 = P2P).
  type ibuf
    integer :: n = 0
    integer, pointer :: t(:), s(:), k(:)
  end type ibuf

contains

  integer function tree_ncells()
    tree_ncells = ncell
  end function tree_ncells

  integer function nterms(p)
    integer, intent(in) :: p
    nterms = (p + 1) * (p + 2) * (p + 3) / 6
  end function nterms

  ! ------------------------------------------------------------- variant

  ! D = 2 always uses the planar operators; D = 3 uses the compressed ones
  ! if asked, else the full ones.
  subroutine select_variant(d, compressed)
    integer, intent(in) :: d
    logical, intent(in) :: compressed
    dim = d
    nch = 2**d
    if (d == 2) then
      variant = V_PLANAR
    else if (compressed) then
      variant = V_COMPRESSED
    else
      variant = V_FULL
    end if
  end subroutine select_variant

  function variant_name()
    character(len=19) :: variant_name
    select case (variant)
    case (V_COMPRESSED)
      variant_name = 'harmonic-compressed'
    case (V_PLANAR)
      variant_name = 'planar'
    case default
      variant_name = 'uncompressed'
    end select
  end function variant_name

  subroutine v_s2m(x, y, z, S, M)
    real(wp), intent(in) :: x, y, z, S(*)
    real(wp), intent(inout) :: M(*)
    select case (variant)
    case (V_COMPRESSED)
      call S2Mc(x, y, z, S, M, order)
    case (V_PLANAR)
      call S2Mxy(x, y, z, S, M, order)
    case default
      call S2M(x, y, z, S, M, order)
    end select
  end subroutine v_s2m

  subroutine v_m2m(x, y, z, M, Ms)
    real(wp), intent(in) :: x, y, z, M(*)
    real(wp), intent(inout) :: Ms(*)
    select case (variant)
    case (V_COMPRESSED)
      call M2Mc(x, y, z, M, Ms, order)
    case (V_PLANAR)
      call M2Mxy(x, y, z, M, Ms, order)
    case default
      call M2M(x, y, z, M, Ms, order)
    end select
  end subroutine v_m2m

  subroutine v_m2l(x, y, z, M, L)
    real(wp), intent(in) :: x, y, z, M(*)
    real(wp), intent(inout) :: L(*)
    select case (variant)
    case (V_COMPRESSED)
      call M2Lc(x, y, z, M, L, order)
    case (V_PLANAR)
      call M2Lxy(x, y, z, M, L, order)
    case default
      call M2L(x, y, z, M, L, order)
    end select
  end subroutine v_m2l

  subroutine v_l2l(x, y, z, L, Ls)
    real(wp), intent(in) :: x, y, z, L(*)
    real(wp), intent(inout) :: Ls(*)
    select case (variant)
    case (V_COMPRESSED)
      call L2Lc(x, y, z, L, Ls, order)
    case (V_PLANAR)
      call L2Lxy(x, y, z, L, Ls, order)
    case default
      call L2L(x, y, z, L, Ls, order)
    end select
  end subroutine v_l2l

  subroutine v_l2p(x, y, z, L, F)
    real(wp), intent(in) :: x, y, z, L(*)
    real(wp), intent(inout) :: F(*)
    select case (variant)
    case (V_COMPRESSED)
      call L2Pc(x, y, z, L, F, order)
    case (V_PLANAR)
      call L2Pxy(x, y, z, L, F, order)
    case default
      call L2P(x, y, z, L, F, order)
    end select
  end subroutine v_l2p

  subroutine v_m2p(x, y, z, M, F)
    real(wp), intent(in) :: x, y, z, M(*)
    real(wp), intent(inout) :: F(*)
    select case (variant)
    case (V_COMPRESSED)
      call M2Pc(x, y, z, M, F, order)
    case (V_PLANAR)
      call M2Pxy(x, y, z, M, F, order)
    case default
      call M2P(x, y, z, M, F, order)
    end select
  end subroutine v_m2p

  ! Sources ib..ie (inclusive) acting on target t, accumulated into the
  ! output of t in Fb. The planar kernel takes no z coordinate.
  subroutine pb(t, ib, ie, Fb)
    integer, intent(in) :: t, ib, ie
    real(wp), intent(inout) :: Fb(*)
    if (ie < ib) return
    if (variant == V_PLANAR) then
      call P2P_batchxy(bx(t), by(t), bx, by, bS, ib, ie, Fb(osz * (t - 1) + 1))
    else
      call P2P_batch(bx(t), by(t), bz(t), bx, by, bz, bS, ib, ie, Fb(osz * (t - 1) + 1))
    end if
  end subroutine pb

  ! Sources in [bo, bo+bn-1] acting on target t, skipping t itself. The
  ! range is split around t so the kernel stays branch-free; excluding by
  ! index (not by r > 0) keeps coincident but distinct particles.
  subroutine p2p_range(t, bo, bn, Fb)
    integer, intent(in) :: t, bo, bn
    real(wp), intent(inout) :: Fb(*)
    if (t >= bo .and. t < bo + bn) then
      call pb(t, bo, t - 1, Fb)
      call pb(t, t + 1, bo + bn - 1, Fb)
    else
      call pb(t, bo, bo + bn - 1, Fb)
    end if
  end subroutine p2p_range

  ! ---------------------------------------------------------------- build

  ! Bit d-1 of the result is set if particle l is above the centre in
  ! dimension d (strictly greater, as in the C++ code).
  integer function octant(l, ctr)
    integer, intent(in) :: l
    real(wp), intent(in) :: ctr(3)
    octant = 0
    if (pos(1, l) > ctr(1)) octant = octant + 1
    if (pos(2, l) > ctr(2)) octant = octant + 2
    if (dim == 3) then
      if (pos(3, l) > ctr(3)) octant = octant + 4
    end if
  end function octant

  ! Split cell c, which holds the particles idx(lo:hi), about its centre
  ! ctr with half-width r. With fill = .false. the cells are only counted,
  ! so the arrays can be sized exactly before the real pass.
  recursive subroutine split(c, lo, hi, ctr, r, fill)
    integer, intent(in) :: c, lo, hi
    real(wp), intent(in) :: ctr(3), r
    logical, intent(in) :: fill
    integer :: i, oc, k, cnt(0:7), start(0:7), nxt(0:7), ch, l
    real(wp) :: cc(3)

    if (fill) then
      nleaf(c) = hi - lo + 1
      boff(c) = lo
    end if
    if (hi - lo + 1 < ncrit) return

    cnt = 0
    do i = lo, hi
      k = octant(idx(i), ctr)
      cnt(k) = cnt(k) + 1
    end do
    start(0) = lo
    do oc = 1, 7
      start(oc) = start(oc - 1) + cnt(oc - 1)
    end do
    nxt = start
    do i = lo, hi
      l = idx(i)
      k = octant(l, ctr)
      tmp(nxt(k)) = l
      nxt(k) = nxt(k) + 1
    end do
    idx(lo:hi) = tmp(lo:hi)

    do oc = 0, nch - 1
      if (cnt(oc) == 0) cycle
      ncell = ncell + 1
      ch = ncell
      cc = ctr
      do k = 1, dim
        cc(k) = ctr(k) + 0.5d0 * r * real(2 * iand(ishft(oc, 1 - k), 1) - 1, wp)
      end do
      if (fill) then
        cen(:, ch) = cc
        rad(ch) = 0.5d0 * r
        rmx(ch) = sqrt(real(dim, wp)) * 0.5d0 * r
        parent(ch) = c
        lvl(ch) = lvl(c) + 1
        child(oc + 1, c) = ch
        nchild(c) = ior(nchild(c), ishft(1, oc))
      end if
      call split(ch, start(oc), start(oc) + cnt(oc) - 1, cc, 0.5d0 * r, fill)
    end do
  end subroutine split

  subroutine buf_push(b, a, s, k)
    type(ibuf), intent(inout) :: b
    integer, intent(in) :: a, s, k
    integer, pointer :: t1(:), t2(:), t3(:)
    integer :: cap
    if (.not. associated(b%t)) then
      allocate(b%t(256), b%s(256), b%k(256))
      b%n = 0
    end if
    cap = size(b%t)
    if (b%n == cap) then
      allocate(t1(2 * cap), t2(2 * cap), t3(2 * cap))
      t1(1:cap) = b%t; t2(1:cap) = b%s; t3(1:cap) = b%k
      deallocate(b%t, b%s, b%k)
      b%t => t1; b%s => t2; b%k => t3
    end if
    b%n = b%n + 1
    b%t(b%n) = a; b%s(b%n) = s; b%k(b%n) = k
  end subroutine buf_push

  ! One step of the dual-tree traversal. Returns 1 = M2L, 2 = P2P,
  ! 3 = expand a, 4 = expand b.
  integer function classify(a, b)
    integer, intent(in) :: a, b
    real(wp) :: dx, dy, dz, r2, rsum
    dx = cen(1, a) - cen(1, b)
    dy = cen(2, a) - cen(2, b)
    dz = cen(3, a) - cen(3, b)
    r2 = dx * dx + dy * dy + dz * dz
    rsum = rmx(a) + rmx(b)
    if (r2 * theta * theta > rsum * rsum) then
      classify = 1
    else if (nchild(a) == 0 .and. nchild(b) == 0) then
      classify = 2
    else if (nchild(b) == 0 .or. (rmx(a) >= rmx(b) .and. nchild(a) /= 0)) then
      classify = 3
    else
      classify = 4
    end if
  end function classify

  recursive subroutine dfs(a, b, buf)
    integer, intent(in) :: a, b
    type(ibuf), intent(inout) :: buf
    integer :: code, oc
    code = classify(a, b)
    select case (code)
    case (1, 2)
      call buf_push(buf, a, b, code)
    case (3)
      do oc = 1, nch
        if (child(oc, a) /= 0) call dfs(child(oc, a), b, buf)
      end do
    case default
      do oc = 1, nch
        if (child(oc, b) /= 0) call dfs(a, child(oc, b), buf)
      end do
    end select
  end subroutine dfs

  ! Dual-tree traversal, in parallel but deterministically: expand
  ! breadth-first until there are enough independent subtree pairs, run each
  ! depth-first into its own buffer, and concatenate the buffers in index
  ! order. The work split and the order are fixed, so the lists are the same
  ! whatever the thread count.
  subroutine traverse()
    type(ibuf) :: head
    type(ibuf), allocatable :: bufs(:)
    integer, allocatable :: fa(:), fb(:), na(:), nb(:)
    integer :: nf, nn, f, oc, code, target_tasks, nthreads, total, pos_, i

    nthreads = 1
!$  nthreads = omp_get_max_threads()
    target_tasks = 8 * nthreads

    nullify(head%t, head%s, head%k)
    allocate(fa(1), fb(1))
    fa(1) = 1; fb(1) = 1
    nf = 1
    do while (nf < target_tasks)
      allocate(na(nch * nf), nb(nch * nf))
      nn = 0
      do f = 1, nf
        code = classify(fa(f), fb(f))
        select case (code)
        case (1, 2)
          call buf_push(head, fa(f), fb(f), code)
        case (3)
          do oc = 1, nch
            if (child(oc, fa(f)) /= 0) then
              nn = nn + 1; na(nn) = child(oc, fa(f)); nb(nn) = fb(f)
            end if
          end do
        case default
          do oc = 1, nch
            if (child(oc, fb(f)) /= 0) then
              nn = nn + 1; na(nn) = fa(f); nb(nn) = child(oc, fb(f))
            end if
          end do
        end select
      end do
      deallocate(fa, fb)
      allocate(fa(max(nn, 1)), fb(max(nn, 1)))
      fa(1:nn) = na(1:nn); fb(1:nn) = nb(1:nn)
      deallocate(na, nb)
      nf = nn
      if (nf == 0) exit
    end do

    allocate(bufs(max(nf, 1)))
    do f = 1, size(bufs)
      nullify(bufs(f)%t, bufs(f)%s, bufs(f)%k)
    end do
!$omp parallel do schedule(dynamic, 1)
    do f = 1, nf
      call dfs(fa(f), fb(f), bufs(f))
    end do
!$omp end parallel do

    total = head%n
    do f = 1, nf
      total = total + bufs(f)%n
    end do
    allocate(il_t(max(total, 1)), il_s(max(total, 1)), il_k(max(total, 1)))
    pos_ = 0
    if (head%n > 0) then
      il_t(1:head%n) = head%t(1:head%n); il_s(1:head%n) = head%s(1:head%n)
      il_k(1:head%n) = head%k(1:head%n)
      pos_ = head%n
    end if
    do f = 1, nf
      i = bufs(f)%n
      if (i > 0) then
        il_t(pos_ + 1:pos_ + i) = bufs(f)%t(1:i)
        il_s(pos_ + 1:pos_ + i) = bufs(f)%s(1:i)
        il_k(pos_ + 1:pos_ + i) = bufs(f)%k(1:i)
        pos_ = pos_ + i
      end if
      if (associated(bufs(f)%t)) deallocate(bufs(f)%t, bufs(f)%s, bufs(f)%k)
    end do
    if (associated(head%t)) deallocate(head%t, head%s, head%k)
    deallocate(bufs, fa, fb)

    call group(1, total, m2l_start, m2l_src)
    call group(2, total, p2p_start, p2p_src)
    deallocate(il_t, il_s, il_k)
  end subroutine traverse

  ! Stable counting sort of the entries of one kind by target cell.
  subroutine group(kind, total, start, src)
    integer, intent(in) :: kind, total
    integer, allocatable, intent(out) :: start(:), src(:)
    integer :: i, c, cnt
    integer, allocatable :: cursor(:)

    allocate(start(ncell + 1), cursor(ncell + 1))
    start = 0
    cnt = 0
    do i = 1, total
      if (il_k(i) == kind) then
        start(il_t(i) + 1) = start(il_t(i) + 1) + 1
        cnt = cnt + 1
      end if
    end do
    start(1) = 1
    do c = 1, ncell
      start(c + 1) = start(c + 1) + start(c)
    end do
    allocate(src(max(cnt, 1)))
    cursor = start
    do i = 1, total
      if (il_k(i) == kind) then
        src(cursor(il_t(i))) = il_s(i)
        cursor(il_t(i)) = cursor(il_t(i)) + 1
      end if
    end do
    deallocate(cursor)
  end subroutine group

  ! position is dim x np (dim from select_variant); strength is ssz x np.
  subroutine build_tree(position, strength, np, ncrit_in, order_in, theta_in)
    integer, intent(in) :: np, ncrit_in, order_in
    real(wp), intent(in) :: position(:,:), strength(ssz, np)
    real(wp), intent(in) :: theta_in
    real(wp) :: avg(3), ext(3), r
    integer :: i, k, maxcells, l, c
    integer, allocatable :: cursor(:)

    n = np; ncrit = ncrit_in; order = order_in; theta = theta_in

    if (allocated(cen)) deallocate(cen, rad, rmx, nleaf, nchild, child, boff, parent, lvl, lev_start, lev_cells)
    if (allocated(bx)) deallocate(bx, by, bz, bS, Fm, perm, Mc, Lc)
    if (allocated(m2l_start)) deallocate(m2l_start, m2l_src, p2p_start, p2p_src)
    if (allocated(idx)) deallocate(idx, tmp, pos)

    allocate(idx(n), tmp(n), pos(3, n))
    pos = 0.0d0
    pos(1:dim, :) = position(1:dim, :)
    do i = 1, n
      idx(i) = i
    end do

    avg = 0.0d0
    do i = 1, n
      avg = avg + pos(:, i)
    end do
    avg = avg / real(n, wp)
    ext = 0.0d0
    do i = 1, n
      do k = 1, dim
        ext(k) = max(ext(k), abs(pos(k, i) - avg(k)))
      end do
    end do
    ! 1.001 so that the root is slightly bigger than the furthest particle.
    r = maxval(ext) * 1.001d0

    ! Pass 1 counts the cells so the arrays can be sized exactly; the
    ! partition it leaves behind is already in tree order, and pass 2
    ! reproduces it.
    ncell = 1
    call split(1, 1, n, avg, r, .false.)
    maxcells = ncell

    allocate(cen(3, maxcells), rad(maxcells), rmx(maxcells), nleaf(maxcells), nchild(maxcells), &
             child(8, maxcells), boff(maxcells), parent(maxcells), lvl(maxcells))
    child = 0; nchild = 0; parent = 0; lvl = 0
    cen(:, 1) = avg
    rad(1) = r
    rmx(1) = sqrt(real(dim, wp)) * r
    ncell = 1
    call split(1, 1, n, avg, r, .true.)

    allocate(bx(n), by(n), bz(n), bS(ssz * n), Fm(osz * n), perm(n))
    bz = 0.0d0
    do i = 1, n
      perm(i) = idx(i)
      bx(i) = pos(1, idx(i))
      by(i) = pos(2, idx(i))
      if (dim == 3) bz(i) = pos(3, idx(i))
      bS(ssz * (i - 1) + 1:ssz * i) = strength(:, idx(i))
    end do

    call traverse()

    ! Cells by level, for the level-synchronous M2M and L2L.
    maxlev = maxval(lvl)
    allocate(lev_start(maxlev + 2), lev_cells(ncell))
    lev_start = 0
    do c = 1, ncell
      lev_start(lvl(c) + 2) = lev_start(lvl(c) + 2) + 1
    end do
    lev_start(1) = 1
    do l = 1, maxlev + 1
      lev_start(l + 1) = lev_start(l + 1) + lev_start(l)
    end do
    allocate(cursor(maxlev + 1))
    cursor = lev_start(1:maxlev + 1)
    do c = 1, ncell
      lev_cells(cursor(lvl(c) + 1)) = c
      cursor(lvl(c) + 1) = cursor(lvl(c) + 1) + 1
    end do
    deallocate(cursor)

    ! Strides come from the variant: the compressed and planar operators use
    ! smaller arrays than Nterms(order).
    select case (variant)
    case (V_COMPRESSED)
      ms = FMMGEN_MULTIPOLESIZE(order); ls = FMMGEN_LOCALSIZE(order)
    case (V_PLANAR)
      ms = FMMGEN_PLANAR_MULTIPOLESIZE(order); ls = FMMGEN_PLANAR_LOCALSIZE(order)
    case default
      ms = nterms(order) - nterms(FMMGEN_SOURCEORDER - 1)
      ls = nterms(order - FMMGEN_SOURCEORDER)
    end select
    allocate(Mc(ms, ncell), Lc(ls, ncell))
    deallocate(idx, tmp, pos)
  end subroutine build_tree

  ! Re-read the source strengths (same ordering as at build time). Positions
  ! and topology never change, so a caller that integrates in time can build
  ! once and call this before each solve.
  subroutine refresh_sources(strength)
    real(wp), intent(in) :: strength(ssz, n)
    integer :: m
!$omp parallel do schedule(static)
    do m = 1, n
      bS(ssz * (m - 1) + 1:ssz * m) = strength(:, perm(m))
    end do
!$omp end parallel do
  end subroutine refresh_sources

  ! ---------------------------------------------------------------- solve

  ! A pure copy through a permutation: race-free, however scheduled.
  subroutine unsort(Fmorton, F)
    real(wp), intent(in) :: Fmorton(*)
    real(wp), intent(out) :: F(osz, n)
    integer :: m
!$omp parallel do schedule(static)
    do m = 1, n
      F(:, perm(m)) = Fmorton(osz * (m - 1) + 1:osz * m)
    end do
!$omp end parallel do
  end subroutine unsort

  subroutine clear_all(clearL)
    logical, intent(in) :: clearL
    integer :: c, m
!$omp parallel do schedule(static)
    do c = 1, ncell
      Mc(:, c) = 0.0d0
    end do
!$omp end parallel do
    if (clearL) then
!$omp parallel do schedule(static)
      do c = 1, ncell
        Lc(:, c) = 0.0d0
      end do
!$omp end parallel do
    end if
!$omp parallel do schedule(static)
    do m = 1, osz * n
      Fm(m) = 0.0d0
    end do
!$omp end parallel do
  end subroutine clear_all

  subroutine p2m_m2m()
    integer :: c, m, l, i, p, oc, ch
!$omp parallel do schedule(static)
    do c = 1, ncell
      if (nleaf(c) < ncrit) then
        do m = boff(c), boff(c) + nleaf(c) - 1
          call v_s2m(cen(1, c) - bx(m), cen(2, c) - by(m), cen(3, c) - bz(m), &
                     bS(ssz * (m - 1) + 1), Mc(1, c))
        end do
      end if
    end do
!$omp end parallel do

    ! Bottom-up, one level at a time, parallel over the PARENTS of the level:
    ! siblings share a parent, so parallelising over children would have
    ! several threads accumulate into the same multipole.
    do l = maxlev, 1, -1
!$omp parallel do schedule(static) private(p, oc, ch)
      do i = lev_start(l), lev_start(l + 1) - 1
        p = lev_cells(i)
        do oc = 1, nch
          ch = child(oc, p)
          if (ch /= 0) then
            call v_m2m(cen(1, p) - cen(1, ch), cen(2, p) - cen(2, ch), cen(3, p) - cen(3, ch), &
                       Mc(1, ch), Mc(1, p))
          end if
        end do
      end do
!$omp end parallel do
    end do
  end subroutine p2m_m2m

  subroutine solve_fmm(F)
    real(wp), intent(out) :: F(osz, n)
    integer :: a, b, i, c, p, m, t, l, oc, ch

    call clear_all(.true.)
    call p2m_m2m()

    ! M2L and P2P are grouped by target, so one thread owns a whole target:
    ! it is the only writer of its L and of its particles' F. No atomics.
!$omp parallel do schedule(dynamic, 16) private(i, b, t)
    do a = 1, ncell
      do i = m2l_start(a), m2l_start(a + 1) - 1
        b = m2l_src(i)
        call v_m2l(cen(1, a) - cen(1, b), cen(2, a) - cen(2, b), cen(3, a) - cen(3, b), &
                   Mc(1, b), Lc(1, a))
      end do
      do i = p2p_start(a), p2p_start(a + 1) - 1
        b = p2p_src(i)
        do t = boff(a), boff(a) + nleaf(a) - 1
          call p2p_range(t, boff(b), nleaf(b), Fm)
        end do
      end do
    end do
!$omp end parallel do

    ! Top-down, one level at a time, parallel over parents (sole writers of
    ! their children's L).
    do l = 1, maxlev
!$omp parallel do schedule(static) private(p, oc, ch)
      do i = lev_start(l), lev_start(l + 1) - 1
        p = lev_cells(i)
        do oc = 1, nch
          ch = child(oc, p)
          if (ch /= 0) then
            call v_l2l(cen(1, ch) - cen(1, p), cen(2, ch) - cen(2, p), cen(3, ch) - cen(3, p), &
                       Lc(1, p), Lc(1, ch))
          end if
        end do
      end do
!$omp end parallel do
    end do

!$omp parallel do schedule(dynamic, 16) private(m)
    do c = 1, ncell
      if (nleaf(c) < ncrit) then
        do m = boff(c), boff(c) + nleaf(c) - 1
          call v_l2p(bx(m) - cen(1, c), by(m) - cen(2, c), bz(m) - cen(3, c), &
                     Lc(1, c), Fm(osz * (m - 1) + 1))
        end do
      end if
    end do
!$omp end parallel do

    call unsort(Fm, F)
  end subroutine solve_fmm

  recursive subroutine bh_walk(m, p)
    integer, intent(in) :: m, p
    integer :: oc, c
    real(wp) :: dx, dy, dz, r2

    if (nleaf(p) >= ncrit) then
      do oc = 1, nch
        c = child(oc, p)
        if (c == 0) cycle
        dx = bx(m) - cen(1, c)
        dy = by(m) - cen(2, c)
        dz = bz(m) - cen(3, c)
        r2 = dx * dx + dy * dy + dz * dz
        ! Squared form of  rad(c) > theta * r.
        if (rad(c) * rad(c) > theta * theta * r2) then
          call bh_walk(m, c)
        else
          call v_m2p(dx, dy, dz, Mc(1, c), Fm(osz * (m - 1) + 1))
        end if
      end do
    else
      call p2p_range(m, boff(p), nleaf(p), Fm)
    end if
  end subroutine bh_walk

  subroutine solve_bh(F)
    real(wp), intent(out) :: F(osz, n)
    integer :: m
    call clear_all(.false.)
    call p2m_m2m()
!$omp parallel do schedule(dynamic, 16)
    do m = 1, n
      call bh_walk(m, 1)
    end do
!$omp end parallel do
    call unsort(Fm, F)
  end subroutine solve_bh

  ! The reference solution. Uses its own scratch buffer, not the persistent
  ! Fm, so it shares no state with the approximate solvers.
  subroutine solve_direct(F)
    real(wp), intent(out) :: F(osz, n)
    real(wp), allocatable :: Fd(:)
    integer :: i
    allocate(Fd(osz * n))
    Fd = 0.0d0
!$omp parallel do schedule(static)
    do i = 1, n
      call p2p_range(i, 1, n, Fd)
    end do
!$omp end parallel do
    call unsort(Fd, F)
    deallocate(Fd)
  end subroutine solve_direct

end module fmm_tree


program fmm_example
  use operators
  use fmm_tree
  implicit none

  integer, parameter :: osz = FMMGEN_OUTPUTSIZE
  integer, parameter :: ssz = FMMGEN_SOURCESIZE
  integer, parameter :: i8 = selected_int_kind(18)
  integer, parameter :: ufld = 21, uerr = 22, utim = 23, upar = 24

  integer :: nparticles = 10000, ncrit = 64, typ = 0, dim = 3, compressed = 0
  integer :: iarg, eq, order, i, k, j
  real(wp) :: theta = 0.4d0
  logical :: nodirect = .false., have_label = .false.
  character(len=256) :: arg, key, val, label = ''
  character(len=512) :: infile = ''
  logical :: have_input = .false.
  character(len=1024) :: line
  integer :: ios
  character(len=512) :: base, fname
  character(len=64) :: fmt
  real(wp), allocatable :: pos(:,:), S(:,:), F_exact(:,:), F_approx(:,:)
  real(wp) :: t_direct, t_approx, errs(osz), l2n(osz), l2d(osz), e, u
  integer, allocatable :: seed(:)
  integer :: nseed

  ! ------------------------------------------------------- command line
  iarg = 1
  do while (iarg <= command_argument_count())
    call get_command_argument(iarg, arg)
    eq = index(arg, '=')
    if (eq > 0) then
      key = arg(1:eq - 1)
      val = arg(eq + 1:)
    else
      key = arg
      val = ''
    end if
    select case (trim(key))
    case ('-d', '--nodirect')
      nodirect = .true.
    case ('-n', '--nparticles', '--ncrit', '-t', '--theta', '--type', &
          '--label', '--compress', '--dim', '--input')
      if (eq == 0) then
        iarg = iarg + 1
        call get_command_argument(iarg, val)
      end if
      select case (trim(key))
      case ('-n', '--nparticles')
        read (val, *) nparticles
      case ('--ncrit')
        read (val, *) ncrit
      case ('-t', '--theta')
        read (val, *) theta
      case ('--type')
        read (val, *) typ
      case ('--label')
        label = val
        have_label = .true.
      case ('--input')
        infile = val
        have_input = .true.
      case ('--compress')
        read (val, *) compressed
      case ('--dim')
        read (val, *) dim
      end select
    case ('-h', '--help')
      print '(a)', 'usage: main [options]'
      print '(a)', '  -n, --nparticles N   total number of particles (10000)'
      print '(a)', '      --ncrit N        maximum number of particles in a cell (64)'
      print '(a)', '  -t, --theta X        opening angle, controls error (0.4)'
      print '(a)', '      --type 0|1       0 = Fast Multipole, 1 = Barnes-Hut (0)'
      print '(a)', '      --label S        label for the output files'
      print '(a)', '      --compress 0|1   harmonic-compressed operators (D=3 only)'
      print '(a)', '      --dim 2|3        2 = planar (z=0, planar operators), 3 = default'
      print '(a)', '      --input FILE     read particles from FILE (one per line: positions then'
      print '(a)', '                       strengths, comma separated) instead of generating them'
      print '(a)', '  -d, --nodirect       skip the direct calculation and the error'
      stop
    case default
      print '(2a)', 'unknown argument: ', trim(arg)
      stop 1
    end select
    iarg = iarg + 1
  end do
  if (typ < 0 .or. typ > 1) stop 'Type must be either 0 (Fast Multipole) or 1 (Barnes-Hut)'
  if (dim /= 2 .and. dim /= 3) stop '--dim must be 2 or 3'

  call select_variant(dim, compressed /= 0)

  ! With --input the particle count comes from the file (one per non-empty
  ! line), overriding --nparticles.
  if (have_input) then
    open (upar, file=trim(infile), status='old', iostat=ios)
    if (ios /= 0) stop 'cannot open input file'
    nparticles = 0
    do
      read (upar, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) nparticles = nparticles + 1
    end do
    close (upar)
  end if

  print '(a)', 'Scaling Test Parameters'
  print '(a)', '-----------------------'
  if (dim == 2) then
    print '(a)', 'Dimension  = 2 (planar)'
  else
    print '(a)', 'Dimension  = 3'
  end if
  print '(a,i0)', 'Nparticles = ', nparticles
  if (have_input) print '(2a)', 'input      = ', trim(infile)
  print '(2a)', 'operators  = ', trim(variant_name())
  print '(a,i0)', 'ncrit      = ', ncrit
  print '(a,f8.6)', 'theta      = ', theta
  print '(a,i0)', 'FMMGEN_MINORDER = ', FMMGEN_MINORDER
  print '(a,i0)', 'FMMGEN_MAXORDER = ', FMMGEN_MAXORDER
  print '(a,i0)', 'FMMGEN_SOURCEORDER = ', FMMGEN_SOURCEORDER
  print '(a,i0)', 'FMMGEN_OUTPUTSIZE = ', FMMGEN_OUTPUTSIZE
  print '(a,i0)', 'FMMGEN_SOURCESIZE = ', FMMGEN_SOURCESIZE
  if (typ == 0) print '(a)', 'FMMGEN TYPE = Fast Multipole Method (Lazy Evaluation)'
  if (typ == 1) print '(a)', 'FMMGEN TYPE = Barnes-Hut Method'

  ! ------------------------------------------------------------ particles
  ! Fixed seed for repeatable runs. random_number is not the C++ example's
  ! generator, so the particles (and hence the errors) differ slightly.
  call random_seed(size=nseed)
  allocate(seed(nseed))
  seed = 12345
  call random_seed(put=seed)

  allocate(pos(dim, nparticles), S(ssz, nparticles))
  allocate(F_exact(osz, nparticles), F_approx(osz, nparticles))
  F_exact = 0.0d0

  if (have_input) then
    ! List-directed reads treat the commas as separators, and the trailing
    ! comma on each line is harmless. The input is not rewritten, so it is
    ! safe to point at particles_n_<N>.txt.
    open (upar, file=trim(infile), status='old')
    i = 0
    do while (i < nparticles)
      read (upar, '(a)') line
      if (len_trim(line) == 0) cycle
      i = i + 1
      read (line, *) pos(:, i), S(:, i)
    end do
    close (upar)
  else
    write (fname, '(a,i0,a)') 'particles_n_', nparticles, '.txt'
    open (upar, file=trim(fname), status='replace')
    write (fmt, '(a,i0,a)') '(', dim + ssz, '(es23.15e3,'',''))'
    do i = 1, nparticles
      do j = 1, dim
        call random_number(u)
        pos(j, i) = (2.0d0 * u - 1.0d0) * 1.0d-9
      end do
      do j = 1, ssz
        call random_number(u)
        S(j, i) = 2.0d0 * u - 1.0d0
      end do
      write (upar, fmt) pos(:, i), S(:, i)
    end do
    close (upar)
  end if

  write (base, '(a,i0,a,i0,a,f8.6,a,i0,a,i0)') '_n_', nparticles, '_ncrit_', ncrit, &
        '_theta_', theta, '_type_', typ, '_d_', dim
  if (have_label) base = trim(base) // '_label_' // trim(label)
  base = trim(base) // '.txt'

  open (utim, file='times' // trim(base), status='replace')

  do order = FMMGEN_MINORDER, FMMGEN_MAXORDER - 1
    call build_tree(pos, S, nparticles, ncrit, order, theta)

    if (order == FMMGEN_MINORDER .and. .not. nodirect) then
      print '(a)', 'Direct'
      print '(a)', '-------'
      t_direct = wtime()
      call solve_direct(F_exact)
      t_direct = wtime() - t_direct
      print '(a,es15.7)', 't_direct = ', t_direct
      write (utim, '(a,es23.15e3)') 'direct,', t_direct
    end if

    print '(a,i0)', 'Order ', order
    print '(a)', '-------'

    F_approx = 0.0d0
    t_approx = wtime()
    if (typ == 0) then
      call solve_fmm(F_approx)
    else
      call solve_bh(F_approx)
    end if
    t_approx = wtime() - t_approx
    write (utim, '(i0,'','',es23.15e3)') order, t_approx

    if (.not. nodirect) then
      write (fmt, '(a,i0,a)') '(', osz, '(es23.15e3,'',''))'
      write (fname, '(a,i0,2a)') 'errors_p_', order, '_', trim(base)
      open (uerr, file=trim(fname), status='replace')
      errs = 0.0d0; l2n = 0.0d0; l2d = 0.0d0
      do i = 1, nparticles
        do k = 1, osz
          e = (F_exact(k, i) - F_approx(k, i)) / F_exact(k, i)
          errs(k) = errs(k) + abs(e)
          l2n(k) = l2n(k) + (F_exact(k, i) - F_approx(k, i))**2
          l2d(k) = l2d(k) + F_exact(k, i)**2
          write (uerr, '(es23.15e3,'','')', advance='no') e
        end do
        write (uerr, '(a)') ''
      end do
      close (uerr)

      ! Same metric as the C++ example (mean |relative error| per component),
      ! plus the L2 relative error, which is not dominated by the few
      ! particles whose exact field happens to be near zero.
      write (fmt, '(a,i0,a)') '(a,', osz, '(es16.7e3,'',''))'
      print fmt, 'Rel errs = ', errs / real(nparticles, wp)
      print fmt, 'L2 errs  = ', sqrt(l2n / l2d)

      write (fname, '(a,i0,2a)') 'field_p_', order, '_', trim(base)
      open (ufld, file=trim(fname), status='replace')
      write (fmt, '(a,i0,a)') '(', 2 * osz, '(es25.16e3,'',''))'
      do i = 1, nparticles
        write (ufld, fmt) (F_exact(k, i), F_approx(k, i), k = 1, osz)
      end do
      close (ufld)
    end if
    print '(a,es15.7,a)', 'Approx. calculation  = ', t_approx, ' seconds.'
  end do
  close (utim)

contains

  ! Wall-clock seconds.
  real(wp) function wtime()
    integer(i8) :: count, rate
    call system_clock(count, rate)
    wtime = real(count, wp) / real(rate, wp)
  end function wtime

end program fmm_example
