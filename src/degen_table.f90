m4_include(util/macro.f90.inc)

module degen_table_m

  use blas95,          only: gemm, dot
  use ieee_arithmetic!, only: ieee_next_after, ieee_positive_inf, ieee_negative_inf
  use lapack95,        only: stev, getrf, getrs, gerfs
  use mpfr_m,          only: mpfr, mpfr_startup, mpfr_cleanup, mpfr_cleanup_cache, add, sub, mul, sqr, fma, fms, fmma, fmms, div, neg, pow, log1p_mpfr => log1p, sqrt_mpfr
  use omp_lib,         only: omp_get_thread_num, omp_get_num_threads
  use zlib_m

  use bin_search_m,    only: bin_search, BS_LESS
  use distributions_m, only: fermi_dirac_integral_1h, &
    &                        fermi_dirac_integral_m1h, &
    &                        fermi_dirac_integral_3h, &
    &                        inv_fermi_dirac_integral_1h
  use error_m,         only: program_error
  use fermi_m,         only: GAMMA, ETA_SMALL, ETA_TINY, ETA_MICRO, ETA_REF, F_REF
  use math_m,          only: linspace, logspace, expm1, log1p, PI
  use qsort_m,         only: qsort
  use util_m,          only: hash
  use vector_m,        only: vector_real, vector_int

  use high_precision_m, hp_log1p => log1p

  implicit none

  integer, parameter :: NGS = 21
  integer, parameter :: MAXN = 2 * NGS - 1
  integer, parameter :: NCHEBY = 5
  integer, parameter :: NCHEBY3 = NCHEBY * NCHEBY * NCHEBY

  real, parameter :: XI_GL(NGS) = [ &
    -0.9937521706203895, -0.9672268385663063, -0.9200993341504008, &
    -0.8533633645833173, -0.7684399634756779, -0.6671388041974123, &
    -0.5516188358872198, -0.4243421202074388, -0.2880213168024011, &
    -0.1455618541608951,  0.0000000000000000,  0.1455618541608951, &
      0.2880213168024011,  0.4243421202074388,  0.5516188358872198, &
      0.6671388041974123,  0.7684399634756779,  0.8533633645833173, &
      0.9200993341504008,  0.9672268385663063,  0.9937521706203895  &
  ]
  real, parameter :: W_GL(NGS) = [ &
    0.0160172282577743, 0.0369537897708525, 0.0571344254268572, &
    0.0761001136283793, 0.0934444234560339, 0.1087972991671484, &
    0.1218314160537285, 0.1322689386333375, 0.1398873947910731, &
    0.1445244039899700, 0.1460811336496904, 0.1445244039899700, &
    0.1398873947910731, 0.1322689386333375, 0.1218314160537285, &
    0.1087972991671484, 0.0934444234560339, 0.0761001136283793, &
    0.0571344254268572, 0.0369537897708525, 0.0160172282577743  &
  ]

  type degen_table
    integer           :: Ndpot  = 0
    integer           :: Ndelta = 0
    integer           :: Neta   = 0
    real              :: eta_max
    real              :: delta_max
    real              :: dpot_max
    real, allocatable :: dpot(:)
    real, allocatable :: delta(:)
    real, allocatable :: eta(:)
    real, allocatable :: pcheby(:,:,:,:)
      !! chebyshev interpolation coefficients per cell
  contains
    procedure :: init       => degen_table_init
    procedure :: load       => degen_table_load
    procedure :: save       => degen_table_save
    procedure :: gen_entry  => degen_table_gen_entry
    procedure :: get        => degen_table_get
    procedure :: get_REG    => degen_table_get_REG
    procedure :: get_INTERP => degen_table_get_INTERP
    procedure :: get_SG     => degen_table_get_SG
    procedure :: get_SGFP   => degen_table_get_SGFP
    procedure :: get_SPLIT1 => degen_table_get_SPLIT1
    procedure :: get_SPLIT2 => degen_table_get_SPLIT2
  end type

  type gauss_table
    real, allocatable :: xi0(:)
    real, allocatable :: xi(:,:)
    real, allocatable :: dxi(:,:)
    real, allocatable :: w(:,:)
    real, allocatable :: dw(:,:)
    real, allocatable :: s0(:)
    real, allocatable :: ds0(:)
  contains
    procedure :: init   => gauss_table_init
    procedure :: gen    => gauss_table_gen
    procedure :: get    => gauss_table_get
    procedure :: interp => gauss_table_interp
    procedure :: output => gauss_table_output
  end type

  type mpfr_vars
    type(mpfr) :: s(0:MAXN), ds(0:MAXN)
    type(mpfr) :: sg(NGS,NGS+1), dsg(NGS,NGS+1)
    type(mpfr) :: xi, ix
    type(mpfr) :: ss, dss
    type(mpfr) :: b, c, dc, d, t, u
    type(mpfr) :: aa(NGS), bb(NGS-1), cc(NGS-1), v(NGS)
    type(mpfr) :: daa(NGS), dbb(NGS-1)
  contains
    procedure :: init     => mpfr_vars_init
    procedure :: destruct => mpfr_vars_destruct
  end type

logical :: DEGEN_TABLE_DEBUG = .false.
public DEGEN_TABLE_DEBUG

contains

  subroutine degen_table_init(this, Meta_lin, Meta_log, eta_lin, eta_max, delta_min, delta_max, Mdpot, dpot_max, folder)
    !! initialize current lookup table for degenerate case
    class(degen_table), intent(out) :: this
    integer,            intent(in)  :: Meta_lin
      !! number of eta cells in table in equidistant region
    integer,            intent(in)  :: Meta_log
      !! number of eta cells in table in logspace region
    real,               intent(in)  :: eta_lin
      !! boundary between eduidistant and logspace region
    real,               intent(in)  :: eta_max
      !! maximal eta on left side of edge
    real,               intent(in)  :: delta_min
      !! minimal delta refinement
    real,               intent(in)  :: delta_max
      !! maximal delta
    integer,            intent(in)  :: Mdpot
      !! number of dpot cells in table
    real,               intent(in)  :: dpot_max
      !! maximal potential drop (minimal drop = 0)
    character(*),       intent(in)  :: folder
      !! folder where to save the table data (without trailing '/')

    integer           :: Neta, Ndelta, Ndpot, counter
    integer           :: ieta, ieta0, ieta1, idelta, idelta0, idelta1, idpot, idpot0, idpot1, jeta, jeta0, jdelta, jdpot
    integer           :: ipiv(NCHEBY3)
    logical           :: status
    character(256)    :: id, fname
    real, allocatable :: eta(:), delta(:), dpot(:), curr(:,:,:)
    real              :: Acheby(NCHEBY3,NCHEBY3), Bcheby(NCHEBY3,NCHEBY3), tmp(NCHEBY3)
    type(gauss_table) :: gtab

    ! try to load
    write (id, "(2I6,4ES25.16E3,I6,ES25.16E3,2I6,3ES25.16E3)") Meta_lin, Meta_log, eta_lin, eta_max, delta_min, delta_max, Mdpot, dpot_max, NGS, NCHEBY, ETA_MICRO, ETA_REF, F_REF
    write (fname, "(2A,Z0,A)") folder, "/degentab_", hash(trim(id)), ".bin"
    call this%load(trim(fname), trim(id), status)
    if (status) return

    ! generate dpot, delta, eta values and chebyshev interpolation
    call gen_axes()

    ! generate gauss lookup table
    call gtab%init(1e-100, 1e5)

    ! memory
    allocate (curr(Neta,-Ndelta:Ndelta,Ndpot), source = 0.0)
    allocate (this%pcheby(NCHEBY3,this%Neta-1,-this%Ndelta:this%Ndelta-1,this%Ndpot-1), source = 0.0)

    counter = 0
    !$omp parallel default(none) &
    !$omp private(ieta,ieta0,ieta1,idelta,idelta0,idelta1,idpot,idpot0,idpot1,jeta,jeta0,jdelta,jdpot,tmp) &
    !$omp shared(this,Neta,Ndelta,Ndpot,counter,ipiv,eta,delta,dpot,curr,Acheby,Bcheby,gtab)

    ! table entries
    !$omp do schedule(dynamic)
    do idpot = 1, Ndpot
      ! delta = 0 => delta_eta = dpot => j = 0
      ! delta > 0
      do idelta = 1, Ndelta
        !$omp atomic
        counter = counter + 1
        !$omp end atomic
        print "(I0,A,I0)", counter, " / ", 2 * Ndpot * Ndelta
        do ieta = 1, Neta
          call this%gen_entry(gtab, dpot(idpot), delta(idelta), eta(ieta), curr(ieta,idelta-1,idpot), curr(ieta,idelta,idpot))
        end do
      end do

      ! delta < 0
      jeta0  = 1
      jdpot = (idpot - 1) / (NCHEBY - 1) + 1
      if (jdpot >= this%Ndpot) jdpot = this%Ndpot - 1
      do idelta = -1, -Ndelta, -1
        !$omp atomic
        counter = counter + 1
        !$omp end atomic
        print "(I0,A,I0,A)", counter, " / ", 2 * Ndpot * Ndelta, " (delta < 0)"

        do while (this%eta(jeta0 + 1) + this%dpot(jdpot + 1) + this%delta((idelta + 1) / (NCHEBY - 1)) < ETA_SMALL)
          jeta0 = jeta0 + 1
        end do
        do ieta = (jeta0 - 1) * (NCHEBY - 1) + 1, Neta
          call this%gen_entry(gtab, dpot(idpot), delta(idelta), eta(ieta), curr(ieta,idelta+1,idpot), curr(ieta,idelta,idpot))
        end do
      end do
    end do
    !$omp end do

    ! interpolation coefficients
    !$omp do schedule(dynamic)
    do jdpot = 1, this%Ndpot - 1
      idpot0 = (jdpot - 1) * (NCHEBY - 1) + 1
      idpot1 = idpot0 + NCHEBY - 1
      do jdelta = 0, this%Ndelta - 1
        idelta0 = jdelta * (NCHEBY - 1)
        idelta1 = idelta0 + NCHEBY - 1
        do jeta = 1, this%Neta - 1
          ieta0 = (jeta - 1) * (NCHEBY - 1) + 1
          ieta1 = ieta0 + NCHEBY - 1

          tmp = reshape(curr(ieta0:ieta1,idelta0:idelta1,idpot0:idpot1), [NCHEBY3])
          this%pcheby(:,jeta,jdelta,jdpot) = tmp
          call getrs(Bcheby, ipiv, this%pcheby(:,jeta,jdelta,jdpot))
          call gerfs(Acheby, Bcheby, ipiv, tmp, this%pcheby(:,jeta,jdelta,jdpot))
        end do
      end do

      jeta0 = 1
      do jdelta = -1, -this%Ndelta, -1
        idelta0 = jdelta * (NCHEBY - 1)
        idelta1 = idelta0 + NCHEBY - 1
        do while (this%eta(jeta0 + 1) + this%dpot(jdpot + 1) + this%delta(jdelta + 1) < ETA_SMALL)
          jeta0 = jeta0 + 1
        end do
        do jeta = jeta0, this%Neta - 1
          ieta0 = (jeta - 1) * (NCHEBY - 1) + 1
          ieta1 = ieta0 + NCHEBY - 1

          tmp = reshape(curr(ieta0:ieta1,idelta0:idelta1,idpot0:idpot1), [NCHEBY3])
          this%pcheby(:,jeta,jdelta,jdpot) = tmp
          call getrs(Bcheby, ipiv, this%pcheby(:,jeta,jdelta,jdpot))
          call gerfs(Acheby, Bcheby, ipiv, tmp, this%pcheby(:,jeta,jdelta,jdpot))
        end do
      end do
    end do
    !$omp end do

    !$omp end parallel

    call this%save(fname, id)

  contains

    subroutine gen_axes()
      integer           :: i, i0, i1, j, k, l, m, n, Mref, Mlin, row, col
      real              :: d, e, t, xcheby(NCHEBY)
      type(vector_real) :: ref

      do i = 1, NCHEBY
        xcheby(i) = cos((PI * (NCHEBY - i)) / (NCHEBY - 1))
      end do

      ! dpot
      this%dpot_max = dpot_max
      this%Ndpot = Mdpot + 1
      this%dpot = linspace(0.0, dpot_max, this%Ndpot)
      Ndpot = Mdpot * (NCHEBY - 1) + 1
      allocate (dpot(Ndpot))
      dpot(1) = this%dpot(1)
      i1 = 1
      do i = 1, Mdpot
        i0 = i1 + 1
        i1 = i1 + NCHEBY - 1
        dpot(i0:i1) = 0.5 * (this%dpot(i) + this%dpot(i+1)) + 0.5 * (this%dpot(i+1) - this%dpot(i)) * xcheby(2:NCHEBY)
      end do

      ! refinement for delta
      call ref%init(0, c = 32)
      call ref%push(0.0)
      d = (this%dpot(Mdpot) - this%dpot(1)) / Mdpot
      e = 1.0
      do while (ref%back() < d)
        t = ref%back() / d
        call ref%push(ref%back() + delta_min * exp(e * t) + (d - delta_min * exp(e)) * expm1(e * t) / expm1(e))
      end do
      Mref = ref%n
      ref%d(1:Mref) = ref%d(1:Mref) * (1 + (d - ref%back()) * ref%d(1:Mref) / ref%back()**2)

      ! delta
      Mlin = nint(delta_max / d) - 1
      this%Ndelta = Mref + Mlin - 1
      allocate (this%delta(-this%Ndelta:this%Ndelta))
      this%delta(0:Mref-1) = ref%d(1:Mref)
      this%delta(Mref:this%Ndelta) = linspace(2 * d, (Mlin + 1) * d, Mlin)
      this%delta_max = this%delta(this%Ndelta)
      this%delta(-this%Ndelta:-1) = -this%delta(this%Ndelta:1:-1)
      Ndelta = this%Ndelta * (NCHEBY - 1)
      allocate (delta(-Ndelta:Ndelta))
      delta(0) = this%delta(0)
      i1 = 0
      do i = 1, this%Ndelta
        i0 = i1 + 1
        i1 = i1 + NCHEBY - 1
        delta(i0:i1) = 0.5 * (this%delta(i-1) + this%delta(i)) + 0.5 * (this%delta(i) - this%delta(i-1)) * xcheby(2:NCHEBY)
        delta(-i0:-i1:-1) = - delta(i0:i1)
      end do

      ! eta
      this%eta_max = eta_max
      this%Neta = Meta_lin + Meta_log + 1
      allocate (this%eta(this%Neta))
      this%eta(1:Meta_lin+1) = linspace(ETA_SMALL, eta_lin, Meta_lin + 1)
      do i = Meta_lin+2, this%Neta
        this%eta(i) = eta_lin * exp(log(eta_max / eta_lin) * (i - 1 - Meta_lin) / (this%Neta - 1 - Meta_lin))
      end do
      Neta = (this%Neta - 1) * (NCHEBY - 1) + 1
      allocate (eta(Neta))
      eta(1) = this%eta(1)
      i1 = 1
      do i = 1, this%Neta - 1
        i0 = i1 + 1
        i1 = i1 + NCHEBY - 1
        eta(i0:i1) = 0.5 * (this%eta(i) + this%eta(i+1)) + 0.5 * (this%eta(i+1) - this%eta(i)) * xcheby(2:NCHEBY)
      end do

      ! chebyshev interpolation
      row = 0
      do k = 1, NCHEBY; do j = 1, NCHEBY; do i = 1, NCHEBY
        row = row + 1
        col = 0
        do n = 0, NCHEBY-1; do m = 0, NCHEBY-1; do l = 0, NCHEBY-1
          col = col + 1
          Acheby(row, col) = xcheby(i)**l * xcheby(j)**m * xcheby(k)**n
        end do; end do; end do
      end do; end do; end do
      Bcheby = Acheby
      call getrf(Bcheby, ipiv = ipiv)
    end subroutine

  end subroutine

  subroutine degen_table_load(this, fname, id, status)
    !! load table from file
    class(degen_table), intent(out) :: this
    character(*),       intent(in)  :: fname
      !! file name
    character(*),       intent(in)  :: id
      !! identifier (arguments passed to init)
    logical,            intent(out) :: status
      !! true => success; false => fail

    integer        :: funit
    character(256) :: id_

    ! check if file exists
    inquire (file = fname, exist = status)
    if (.not. status) return

    open (newunit = funit, file = fname, status = "old", action = "read", form = "unformatted")

    read (funit) id_
    if (id_ /= id) then
      status = .false.
      return
    end if

    read (funit) this%Ndpot
    read (funit) this%Ndelta
    read (funit) this%Neta
    read (funit) this%eta_max
    read (funit) this%delta_max
    read (funit) this%dpot_max
    allocate (this%dpot(this%Ndpot), this%delta(-this%Ndelta:this%Ndelta), this%eta(this%Neta))
    allocate (this%pcheby(NCHEBY3,this%Neta-1,-this%Ndelta:this%Ndelta-1,this%Ndpot-1), source = 0.0)
    read (funit) this%dpot
    read (funit) this%delta
    read (funit) this%eta
    read (funit) this%pcheby

    close (funit)
    status = .true.
  end subroutine

  subroutine degen_table_save(this, fname, id)
    !! save table to file
    class(degen_table), intent(in) :: this
    character(*),       intent(in) :: fname
      !! file name
    character(*),       intent(in) :: id
      !! identifier (arguments passed to init)

    integer :: funit

    open (newunit = funit, file = fname, status = "replace", action = "write", form = "unformatted")
    write (funit) id
    write (funit) this%Ndpot
    write (funit) this%Ndelta
    write (funit) this%Neta
    write (funit) this%eta_max
    write (funit) this%delta_max
    write (funit) this%dpot_max
    write (funit) this%dpot
    write (funit) this%delta
    write (funit) this%eta
    write (funit) this%pcheby
    close (funit)
  end subroutine

  subroutine degen_table_gen_entry(this, gtab, dpot, delta, eta1, curr0, curr)
    !! generate one entry in table using gaussian quadrature
    use math_m, only: ber

    class(degen_table), intent(in)    :: this
    type(gauss_table),  intent(in)    :: gtab
      !! gauss lookup table
    real,               intent(in)    :: dpot
      !! potential drop
    real,               intent(in)    :: delta
      !! delta = eta2 - eta1 - dpot
    real,               intent(in)    :: eta1
      !! eta at left end of edge
    real,               intent(in)    :: curr0
      !! initial guess for edge current (from previous table entry)
    real,               intent(out)   :: curr
      !! output new edge current

    real, parameter :: RTOL = 5e-16

    integer :: it
    logical :: lpole, rpole
    real    :: delta_eta, delta_eta_sgn, eta2, eta_min, eta_max, dens1, dens2, ddens1, ddens2, curr_min, curr_max, curr_old
    real    :: f, df, dcurr, err

    m4_ignore(this)

    delta_eta = dpot + delta
    eta2      = eta1 + delta_eta

    ! special case: dpot == 0 => curr = - int_{eta1}^{eta2} F12(eta) deta
    if (dpot == 0) then
      ! int F12(eta) deta = F32(eta)
      call fermi_dirac_integral_3h(eta1, dens1, ddens1)
      call fermi_dirac_integral_3h(eta2, dens2, ddens2)
      curr = dens1 - dens2
      return
    end if

    ! special case: eta = const => deta/dx = dpot - j / F12(eta) = 0
    if (abs(delta_eta) < 1e-6) then
      call fermi_dirac_integral_1h( eta1, dens1, ddens1)
      call fermi_dirac_integral_m1h(eta1, dens2, ddens2)
      curr = dens1 * (dpot - delta_eta * ber(dpot * dens2 / dens1))
      return
    end if

    ! get densities
    call fermi_dirac_integral_1h(eta1, dens1, ddens1)
    call fermi_dirac_integral_1h(eta2, dens2, ddens2)

    eta_min       = min(eta1, eta2)
    eta_max       = max(eta1, eta2)
    delta_eta_sgn = sign(1.0, delta_eta)
    delta_eta     = abs(delta_eta)

    ! get current range
    call get_range()

    ! remove possible infinities at interval bounds
    if (lpole) then
      call residual(curr_min, f, df)
      do while (.not. ieee_is_finite(f))
        curr_min = ieee_next_after(curr_min, huge(1.0))
        call residual(curr_min, f, df)
      end do
      if (f >= 0) then
        ! solution already found (sign change from -1 to +1)
        curr = curr_min
        return
      end if
    end if
    if (rpole) then
      call residual(curr_max, f, df)
      do while (.not. ieee_is_finite(f))
        curr_max = ieee_next_after(curr_max, -huge(1.0))
        call residual(curr_max, f, df)
      end do
      if (f <= 0) then
        ! solution already found (sign change from +1 to -1)
        curr = curr_max
        return
      end if
    end if

    ! initial guess
    curr = curr0
    if (curr < curr_min) then
      curr = curr_min
    elseif (curr > curr_max) then
      curr = curr_max
    end if

    ! newton iteration / bisection
    err = 0.5 * (curr_max - curr_min)
    it  = 0
    do while (((curr_max - curr_min) > 0.5 * abs(curr) * RTOL) .and. (err > abs(curr) * RTOL))
      it = it + 1

      ! evaluate residual
      call residual(curr, f, df)
      do while (.not. ieee_is_finite(f))
        if (lpole) then
          curr     = ieee_next_after(curr, huge(1.0))
          curr_min = curr
        else
          curr     = ieee_next_after(curr, -huge(1.0))
          curr_max = curr
        end if
        call residual(curr, f, df)
      end do

      ! newton update
      dcurr = - f / df
      err   = abs(dcurr)

      ! update bounds (assume monotonic behaviour)
      if (dcurr > 0) then
        curr_min = curr
      else
        curr_max = curr
      end if

      ! update solution
      curr_old = curr
      curr = curr + dcurr

      ! bisection
      if ((curr < curr_min) .or. (curr > curr_max) .or. ((curr_old == curr_min) .and. (curr == curr_max))) then
        curr = 0.5 * (curr_min + curr_max)
      end if
    end do

  contains

    subroutine get_range()
      real :: tmp, densc

      ! get current min, max by slope (detadx must be equal to eta2 - eta1 for some x in [0, 1])
      curr_min = min(dens1, dens2) * abs(delta)
      curr_max = max(dens1, dens2) * abs(delta)
      if (sign(1.0, delta) > 0) then
        ! correct sign and swap curr_min, curr_max
        tmp      = curr_min
        curr_min = - curr_max
        curr_max = - tmp
      end if

      ! further reduce current range if possible
      lpole = .false.
      rpole = .false.
      if (eta2 < eta1) then
        curr_min = max(curr_min, dpot * dens1)
        lpole = (curr_min == dpot * dens1)
      elseif (eta2 <= eta1 + dpot) then
        curr_min = max(curr_min, 0.0)
        curr_max = min(curr_max, dpot * dens1)
        rpole = (curr_max == dpot * dens1)
      else
        call fermi_dirac_integral_1h(eta2 - 0.5 * dpot, densc, tmp)
        curr_min = max(curr_min, - densc * (delta_eta - dpot))
        curr_max = min(curr_max, 0.0)
      end if
    end subroutine

    subroutine residual(curr, f, df)
      real, intent(in)  :: curr
      real, intent(out) :: f
      real, intent(out) :: df

      integer       :: k
      real          :: xi0, dxi0, eta, eta0, deta, deta0, dens, ddens, xi(NGS), dxi(NGS), w(NGS), dw(NGS), s0, ds0
      type(hp_real) :: hcurr, hf, dhf, hg, dhg, ht, dht, hu, dhu

      hcurr = real_to_hp(curr)
      hf    = real_to_hp(- delta_eta_sgn)
      dhf   = real_to_hp(0.0)

      if (curr / dpot <= 0) then
        ! Gauss-Legendre quadrature
        xi  = 0.5 * (XI_GL + 1.0)
        w   = 0.5 * W_GL
        do k = 1, NGS
          eta = eta_min + delta_eta * xi(k)
          call fermi_dirac_integral_1h(eta, dens, ddens)
          ht = 1.0 / (TwoProduct(dpot, dens) - hcurr)
          hf   =  hf + delta_eta * w(k) * dens * ht
          dhf  = dhf + delta_eta * w(k) * dens * ht**2
        end do
      else
        call inv_fermi_dirac_integral_1h(curr / dpot, eta0, deta0)
        deta0 = deta0 / dpot

        if (eta0 < eta_min) then
          xi0  = (eta_min - eta0) / delta_eta
          dxi0 = - deta0 / delta_eta
        elseif (eta0 > eta_max) then
          xi0  = (eta0 - eta_max) / delta_eta
          dxi0 = deta0 / delta_eta
        else
          ! pole within interval => f = +/- infinity depending on position of pole
          if (lpole) then
            f  = ieee_value(f, ieee_negative_inf)
            df = ieee_value(f, ieee_positive_inf)
          else
            f  = ieee_value(f, ieee_positive_inf)
            df = ieee_value(f, ieee_positive_inf)
          end if
          return
        end if

        call gtab%get(xi0, xi, dxi, w, dw, s0, ds0)
        dxi = dxi * dxi0
        dw  = dw  * dxi0
        ds0 = ds0 * dxi0

        do k = 1, NGS
          if (eta0 < eta_min) then
            eta  = eta_min + delta_eta *  xi(k)
            deta =           delta_eta * dxi(k)
          else
            eta  = eta_max - delta_eta * xi(k)
            deta =         - delta_eta * dxi(k)
          end if

          call fermi_dirac_integral_1h(eta, dens, ddens)
          ddens = ddens * deta

          ht  = 1.0 / (TwoProduct(dpot, dens) - hcurr)
          dht = ht**2 * (1 - TwoProduct(dpot, ddens))

          hu  = TwoSum( eta0, -  eta)
          dhu = TwoSum(deta0, - deta)

          if (eta0 < eta_min) then
            hu  = - hu
            dhu = - dhu
          end if

          hg  = dens * ht * hu
          dhg = ddens * ht * hu + dens * dht * hu + dens * ht * dhu

          hf  = hf + (w(k) * hg) * s0
          dhf = dhf + (dw(k) * hg + w(k) * dhg) * s0 + (w(k) * hg) * ds0
        end do
      end if

      f  = hp_to_real(hf)
      df = hp_to_real(dhf)
    end subroutine

  end subroutine

  subroutine degen_table_get(this, eta, dpot, j, djdeta, djddpot)
    !! lookup current
    class(degen_table), intent(in)  :: this
    real,               intent(in)  :: eta(2)
      !! left/right eta
    real,               intent(in)  :: dpot
      !! normalized potential drop
    real,               intent(out) :: j
      !! output normalized edge current
    real,               intent(out) :: djdeta(2)
      !! output derivatives of j wrt eta
    real,               intent(out) :: djddpot
      !! output derivatives of j wrt dpot

    logical :: flip
    real    :: eta_(2), dpot_

    ! flip edge direction if potential drop is negative
    flip = (dpot < 0)
    if (flip) then
      eta_(1) = eta(2)
      eta_(2) = eta(1)
      dpot_   = - dpot
    else
      eta_  = eta
      dpot_ = dpot
    end if

    if (all(eta <= ETA_MICRO)) then
      ! if (DEGEN_TABLE_DEBUG) print "(A,3ES25.16E3,A)", "REG:    ", eta, dpot
      ! use regularization to avoid numbers that are too small
      call this%get_REG(eta_, dpot_, j, djdeta, djddpot)
    elseif (all(eta >= ETA_SMALL)) then
      ! if (DEGEN_TABLE_DEBUG) print "(A,3ES25.16E3,A)", "INTERP:    ", eta, dpot
      ! directly interpolate from table
      call this%get_INTERP(eta_, dpot_, j, djdeta, djddpot)
    elseif (any(eta > ETA_SMALL)) then
      if (any(eta < ETA_MICRO)) then
        if (DEGEN_TABLE_DEBUG) print "(A,3ES25.16E3,A)", "SPLIT2:    ", eta, dpot
        ! split between regularized (R), approximated (S) and full Fermi-Dirac-Integral (T) regions (RS and ST split)
        call this%get_SPLIT2(eta_, dpot_, j, djdeta, djddpot)
      else
        if (DEGEN_TABLE_DEBUG) print "(A,3ES25.16E3,A)", "SPLIT1 ST:    ", eta, dpot
        ! split between approximated (S) and full Fermi-Dirac-Integral (T) regions
        call this%get_SPLIT1(eta_, dpot_, j, djdeta, djddpot)
      end if
    elseif (any(eta < ETA_MICRO)) then
      if (DEGEN_TABLE_DEBUG) print "(A,3ES25.16E3,A)", "SPLIT1 RS:    ", eta, dpot
      ! split between regularized (R) and approximated Fermi-Dirac-Integral (S) regions
      call this%get_SPLIT1(eta_, dpot_, j, djdeta, djddpot)
    elseif (all(eta <= ETA_TINY)) then
      ! if (DEGEN_TABLE_DEBUG) print "(A,3ES25.16E3,A)", "SG:    ", eta, dpot
      ! Scharfetter-Gummel
      call this%get_SG(eta_, dpot_, j, djdeta, djddpot)
    else
      ! if (DEGEN_TABLE_DEBUG) print "(A,3ES25.16E3,A)", "SGFP:    ", eta, dpot
      ! modified Scharfetter-Gummel fixed-point iteration (approximated Fermi-Dirac-Integral region)
      call this%get_SGFP(eta_, dpot_, j, djdeta, djddpot)
    end if

    ! flip edge direction back (dj/ddpot = d(-j)/d(-dpot) unchanged)
    if (flip) then
      j = - j
      djdeta = - [djdeta(2), djdeta(1)]
    end if
  end subroutine

  subroutine degen_table_get_REG(this, eta, dpot, j, djdeta, djddpot)
    !! use regularization to get edge current
    use math_m, only: ber, dberdx

    class(degen_table), intent(in)  :: this
    real,               intent(in)  :: eta(2)
      !! left/right eta
    real,               intent(in)  :: dpot
      !! normalized potential drop
    real,               intent(out) :: j
      !! output interpolated edge current
    real,               intent(out) :: djdeta(2)
      !! output derivatives of j wrt eta
    real,               intent(out) :: djddpot
      !! output derivatives of j wrt dpot

    real, parameter :: beta  = (ETA_MICRO - log(F_REF)) / (ETA_MICRO - ETA_REF)
    real, parameter :: alpha = exp(beta * ETA_REF) / F_REF

    real :: B1, B2, n1, n2, h

    m4_ignore(this)

    B1 = ber( beta * dpot)
    B2 = ber(-beta * dpot)

    n1 = exp(beta * eta(1)) / alpha
    n2 = exp(beta * eta(2)) / alpha

    h = hp_to_real(TwoSum(eta(2), -eta(1)) - dpot)
    j = - B2 * n1 * expm1(beta * h) / beta

    djdeta(1) =   B2 * n1
    djdeta(2) = - B1 * n2
    djddpot   = - dberdx(-beta * dpot) * n1 - dberdx(beta * dpot) * n2
  end subroutine

  subroutine degen_table_get_INTERP(this, eta, dpot, j, djdeta, djddpot)
    !! interpolate current
    class(degen_table), intent(in)  :: this
    real,               intent(in)  :: eta(2)
      !! left/right eta
    real,               intent(in)  :: dpot
      !! normalized potential drop
    real,               intent(out) :: j
      !! output interpolated edge current
    real,               intent(out) :: djdeta(2)
      !! output derivatives of j wrt eta
    real,               intent(out) :: djddpot
      !! output derivatives of j wrt dpot

    integer :: jeta, jdelta, jdpot, k, l, m, n
    real    :: delta, u, du, v, dv, w, dw, t(NCHEBY3), dtdu(NCHEBY3), dtdv(NCHEBY3), dtdw(NCHEBY3), djdu, djdv, djdw

    delta = eta(2) - eta(1) - dpot

    if ((eta(1) < ETA_SMALL) .or. (eta(1) > this%eta_max)) then
      print "(A,ES25.16E3)", "eta(1) = ", eta(1)
      call program_error("eta(1) out of range")
    elseif (eta(1) == ETA_SMALL) then
      jeta = 1
    elseif (eta(1) == this%eta_max) then
      jeta = this%Neta - 1
    else
      jeta = bin_search(this%eta,   eta(1), BS_LESS)
    end if

    if (abs(delta) > this%delta_max) then
      print "(A,ES25.16E3)", "delta = ", delta
      call program_error("delta out of range")
    elseif (delta == - this%delta_max) then
      jdelta = - this%Ndelta
    elseif (delta == this%delta_max) then
      jdelta = this%Ndelta - 1
    else
      jdelta = bin_search(this%delta, delta,  BS_LESS) - (size(this%delta) - 1) / 2 - 1
    end if

    if (dpot > this%dpot_max) then
      print "(A,ES25.16E3)", "dpot = ", dpot
      call program_error("dpot out of range")
    elseif (dpot == 0) then
      jdpot = 1
    elseif (dpot == this%dpot_max) then
      jdpot = this%Ndpot - 1
    else
      jdpot  = bin_search(this%dpot,  dpot,   BS_LESS)
    end if

    u = (2 * eta(1) - this%eta(  jeta  ) - this%eta(  jeta  +1)) / (this%eta(  jeta  +1) - this%eta(  jeta  ))
    v = (2 * delta  - this%delta(jdelta) - this%delta(jdelta+1)) / (this%delta(jdelta+1) - this%delta(jdelta))
    w = (2 * dpot   - this%dpot( jdpot ) - this%dpot( jdpot +1)) / (this%dpot( jdpot +1) - this%dpot( jdpot ))

    du = 2.0 / (this%eta(  jeta  +1) - this%eta(  jeta  ))
    dv = 2.0 / (this%delta(jdelta+1) - this%delta(jdelta))
    dw = 2.0 / (this%dpot( jdpot +1) - this%dpot( jdpot ))

    n = 0
    do m = 0, NCHEBY-1; do l = 0, NCHEBY-1; do k = 0, NCHEBY-1
      n = n + 1
      t(n) = u**k * v**l * w**m
      if (k > 0) then
        dtdu(n) = k * u**(k-1) * v**l * w**m
      else
        dtdu(n) = 0
      end if
      if (l > 0) then
        dtdv(n) = l * u**k * v**(l-1) * w**m
      else
        dtdv(n) = 0
      end if
      if (m > 0) then
        dtdw(n) = m * u**k * v**l * w**(m-1)
      else
        dtdw(n) = 0
      end if
    end do; end do; end do
    j    = dot(t, this%pcheby(:,jeta,jdelta,jdpot))
    djdu = dot(dtdu, this%pcheby(:,jeta,jdelta,jdpot))
    djdv = dot(dtdv, this%pcheby(:,jeta,jdelta,jdpot))
    djdw = dot(dtdw, this%pcheby(:,jeta,jdelta,jdpot))

    djdeta(1) = djdu * du - djdv * dv
    djdeta(2) = djdv * dv
    djddpot   = djdw * dw - djdv * dv
  end subroutine

  subroutine degen_table_get_SG(this, eta, dpot, j, djdeta, djddpot)
    !! get Scharfetter-Gummel Current
    use math_m, only: ber, dberdx

    class(degen_table), intent(in)  :: this
    real,               intent(in)  :: eta(2)
      !! left/right eta
    real,               intent(in)  :: dpot
      !! normalized potential drop
    real,               intent(out) :: j
      !! output interpolated edge current
    real,               intent(out) :: djdeta(2)
      !! output derivatives of j wrt eta
    real,               intent(out) :: djddpot
      !! output derivatives of j wrt dpot

    real :: B1, B2, n1, n2

    m4_ignore(this)

    B1 = ber( dpot)
    B2 = ber(-dpot)

    n1 = exp(eta(1))
    n2 = exp(eta(2))

    j = - B2 * n1 * hp_to_real(expm1(TwoSum(eta(2), -eta(1)) - dpot))

    djdeta(1) =   B2 * n1
    djdeta(2) = - B1 * n2
    djddpot   = - dberdx(-dpot) * n1 - dberdx(dpot) * n2
  end subroutine

  subroutine degen_table_get_SGFP(this, eta, dpot, j, djdeta, djddpot)
    !! Scharfetter-Gummel fixed-point iteration
    use math_m, only: ber1 => ber, dberdx

    class(degen_table), intent(in)  :: this
    real,               intent(in)  :: eta(2)
      !! left/right eta
    real,               intent(in)  :: dpot
      !! normalized potential drop
    real,               intent(out) :: j
      !! output interpolated edge current
    real,               intent(out) :: djdeta(2)
      !! output derivatives of j wrt eta
    real,               intent(out) :: djddpot
      !! output derivatives of j wrt dpot

    real, parameter :: RTOL = 1e-12

    real          :: n1, n2, j0, b, c
    type(hp_real) :: h

    m4_ignore(this)

    n1 = exp(eta(1))
    n2 = exp(eta(2))

    ! start with pure Scharfetter-Gummel
    h = TwoSum(eta(2), -eta(1)) - dpot
    j = - n1 * ber1(-dpot) * hp_to_real(expm1(h))

    ! fixed point iteration
    j0 = huge(1.0)
    do while (j /= j0)
      j0 = j
      j  = - n1 * hp_to_real(ber(TwoSum(-dpot, GAMMA * j0)) * expm1(h + GAMMA * j0))
    end do

    ! derivatives
    b         = n2 * dberdx( dpot - GAMMA * j) + n1 * dberdx(-dpot + GAMMA * j)
    c         = 1.0 / (1.0 - GAMMA * b)
    djdeta(1) =   n1 * ber1(-dpot + GAMMA * j) * c
    djdeta(2) = - n2 * ber1( dpot - GAMMA * j) * c
    djddpot   = - b * c
  end subroutine

  subroutine degen_table_get_SPLIT1(this, eta, dpot, j, djdeta, djddpot)
    !! split edge at single point to evaluate current
    class(degen_table), intent(in)  :: this
    real,               intent(in)  :: eta(2)
      !! left/right eta
    real,               intent(in)  :: dpot
      !! normalized potential drop
    real,               intent(out) :: j
      !! output interpolated edge current
    real,               intent(out) :: djdeta(2)
      !! output derivatives of j wrt eta
    real,               intent(out) :: djddpot
      !! output derivatives of j wrt dpot

    real, parameter :: RTOL = 1e-14

    integer :: it, split
    real    :: a, eta_split, x, xmin, xmax, xold, dx, dxdeta(2), dxddpot, err, f, dfdx, dfdeta(2), dfddpot

    ! allowed range
    xmin = 0
    xmax = 1

    ! split at which eta?
    if ((eta(1) < ETA_MICRO) .or. (eta(2) < ETA_MICRO)) then
      split     = 1 ! RS
      eta_split = ETA_MICRO
    else
      split     = 2 ! ST
      eta_split = ETA_SMALL
    end if

    ! initial guess using Scharfetter-Gummel formula for eta(x) (not a good approximation for REG-SG split)
    a = (exp(eta_split) - exp(eta(1))) / (exp(eta(2)) - exp(eta(1)))
    if (abs(dpot) < 1e-3) then
      x = a - 0.5 * a * (a - 1) * dpot
    else
      x = log1p(a * expm1(dpot)) / dpot
    end if
    if (x < 0) x = 0
    if (x > 1) x = 1

    it = 0
    err = huge(1.0)
    do while ((it < 2) .or. (((xmax - xmin) > 0.5 * abs(x) * RTOL) .and. (err > abs(x) * RTOL)))
      it = it + 1

      if (it > 25) then
        print "(A,ES25.16E3,A,ES25.16E3,A)", "eta = [", eta(1), ", ", eta(2), "]"
        print "(A,ES25.16E3)", "dpot = ", dpot
        call program_error("No convergence")
      end if

      ! Newton update
      call residual(x, f, dfdx, dfdeta, dfddpot)
      dx  = -f / dfdx
      err = abs(dx)

      ! update bounds (assume monotonic behaviour)
      if (dx > 0) then
        xmin = x
      else
        xmax = x
      end if

      ! print "(A,I0,A,ES25.16E3,A,ES25.16E3)", "  it = ", it, "; x = ", x, "; dx = ", dx

      ! update solution
      xold = x
      x    = x + dx

      ! bisection
      if ((x < xmin) .or. (x > xmax) .or. ((xold == xmin) .and. (x == xmax))) then
        x = 0.5 * (xmin + xmax)
      end if
    end do

    ! implicit differentiation
    dxdeta  = - dfdeta / dfdx
    dxddpot = - dfddpot / dfdx

    ! evaluate current
    call eval_j(x, dxdeta, dxddpot, j, djdeta, djddpot)

  contains

    subroutine residual(x, f, dfdx, dfdeta, dfddpot)
      !! difference between the current in the left and right part of the edge
      real, intent(in)  :: x
        !! edge split position (solution variable)
      real, intent(out) :: f
        !! output residual: j1*(1-x) - j2*x
      real, intent(out) :: dfdx
        !! output derivative of f wrt x
      real, intent(out) :: dfdeta(2)
        !! output derivative of f wrt eta
      real, intent(out) :: dfddpot
        !! output derivative of f wrt dpot

      real :: dpot1, dpot2, j1, j2, dj1deta(2), dj2deta(2), dj1ddpot1, dj2ddpot2

      ! get current in left part
      dpot1 = dpot * x
      if (split == 1) then ! RS
        if (eta(1) < ETA_MICRO) then
          call this%get_REG([eta(1), eta_split], dpot1, j1, dj1deta, dj1ddpot1)
        else
          call this%get_SGFP([eta(1), eta_split], dpot1, j1, dj1deta, dj1ddpot1)
        end if
      else ! ST
        if (eta(1) < ETA_SMALL) then
          call this%get_SGFP([eta(1), eta_split], dpot1, j1, dj1deta, dj1ddpot1)
        else
          call this%get_INTERP([eta(1), eta_split], dpot1, j1, dj1deta, dj1ddpot1)
        end if
      end if

      ! get current in right part
      dpot2 = dpot * (1 - x)
      if (split == 1) then ! RS
        if (eta(2) < ETA_MICRO) then
          call this%get_REG([eta_split, eta(2)], dpot2, j2, dj2deta, dj2ddpot2)
        else
          call this%get_SGFP([eta_split, eta(2)], dpot2, j2, dj2deta, dj2ddpot2)
        end if
      else ! ST
        if (eta(2) < ETA_SMALL) then
          call this%get_SGFP([eta_split, eta(2)], dpot2, j2, dj2deta, dj2ddpot2)
        else
          call this%get_INTERP([eta_split, eta(2)], dpot2, j2, dj2deta, dj2ddpot2)
        end if
      end if

      ! residual (respect scale)
      f         = j1 * (1 - x) - j2 * x ! j1 / x == j2 / (1 - x)
      dfdx      = - j1 - j2 + dj1ddpot1 * dpot2 + dj2ddpot2 * dpot1
      dfdeta(1) =   dj1deta(1) * (1 - x)
      dfdeta(2) = - dj2deta(2) * x
      dfddpot   = (dj1ddpot1 - dj2ddpot2) * x * (1 - x)
    end subroutine

    subroutine eval_j(x, dxdeta, dxddpot, j, djdeta, djddpot)
      !! evaluate current after splitting position has been found
      real, intent(in)  :: x
      real, intent(in)  :: dxdeta(2)
      real, intent(in)  :: dxddpot
      real, intent(out) :: j
      real, intent(out) :: djdeta(2)
      real, intent(out) :: djddpot

      real :: dpot1, j1, dj1deta(2), dj1ddpot1, dj1dx
      real :: dpot2, j2, dj2deta(2), dj2ddpot2, dj2dx

      ! evaluate current using longer part of edge (avoid division by small number)
      if (x >= 0.5) then
        dpot1 = dpot * x
        if (split == 1) then ! RS
          if (eta(1) < ETA_MICRO) then
            call this%get_REG([eta(1), eta_split], dpot1, j1, dj1deta, dj1ddpot1)
          else
            call this%get_SGFP([eta(1), eta_split], dpot1, j1, dj1deta, dj1ddpot1)
          end if
        else ! ST
          if (eta(1) < ETA_SMALL) then
            call this%get_SGFP([eta(1), eta_split], dpot1, j1, dj1deta, dj1ddpot1)
          else
            call this%get_interp([eta(1), eta_split], dpot1, j1, dj1deta, dj1ddpot1)
          end if
        end if

        j         = j1 / x
        dj1dx     = dj1ddpot1 * dpot
        djdeta(1) = (dj1deta(1) + dxdeta(1) * (dj1dx - j)) / x
        djdeta(2) =               dxdeta(2) * (dj1dx - j)  / x
        djddpot   = dj1ddpot1   + dxddpot   * (dj1dx - j)  / x
      else
        dpot2 = dpot * (1 - x)
        if (split == 1) then ! RS
          if (eta(2) < ETA_MICRO) then
            call this%get_REG([eta_split, eta(2)], dpot2, j2, dj2deta, dj2ddpot2)
          else
            call this%get_SGFP([eta_split, eta(2)], dpot2, j2, dj2deta, dj2ddpot2)
          end if
        else ! ST
          if (eta(2) < ETA_SMALL) then
            call this%get_SGFP([eta_split, eta(2)], dpot2, j2, dj2deta, dj2ddpot2)
          else
            call this%get_interp([eta_split, eta(2)], dpot2, j2, dj2deta, dj2ddpot2)
          end if
        end if

        j         = j2 / (1 - x)
        dj2dx     = - dj2ddpot2 * dpot
        djdeta(1) =               dxdeta(1) * (dj2dx + j)  / (1 - x)
        djdeta(2) = (dj2deta(2) + dxdeta(2) * (dj2dx + j)) / (1 - x)
        djddpot   = dj2ddpot2   + dxddpot   * (dj2dx + j)  / (1 - x)
      end if
    end subroutine

  end subroutine

  subroutine degen_table_get_SPLIT2(this, eta, dpot, j, djdeta, djddpot)
    !! split edge at two points to evaluate current
    class(degen_table), intent(in)  :: this
    real,               intent(in)  :: eta(2)
      !! left/right eta
    real,               intent(in)  :: dpot
      !! normalized potential drop
    real,               intent(out) :: j
      !! output interpolated edge current
    real,               intent(out) :: djdeta(2)
      !! output derivatives of j wrt eta
    real,               intent(out) :: djddpot
      !! output derivatives of j wrt dpot

    real, parameter :: RTOL = 1e-12

    integer :: it, i, ipiv(2)
    real    :: a(2), eta_split(2), x(2), xold(2), dx(2), dxdeta(2,2), dxddpot(2), err, err_dx(2), djdx(2)
    real    :: f(2), dfdx(2,2), dfdeta(2,2), dfddpot(2), lu(2,2)

    ! initial guess using Scharfetter-Gummel formula for eta(x)
    if (eta(1) < eta(2)) then
      eta_split = [ETA_MICRO, ETA_SMALL]
    else
      eta_split = [ETA_SMALL, ETA_MICRO]
    end if

    a = (exp(eta_split) - exp(eta(1))) / (exp(eta(2)) - exp(eta(1)))
    ! a = expm1(eta_split - eta(1)) / expm1(eta(2) - eta(1))
    if (abs(dpot) < 1e-3) then
      x = a - 0.5 * a * (a - 1) * dpot
    elseif (abs(dpot) < 500) then
      x = log1p(a * expm1(dpot)) / dpot
    else
      x = 1 + log(a) / dpot
    end if
    do i = 1, 2
      if (x(i) < 0) x(i) = 0
      if (x(i) > 1) x(i) = 1
    end do

    ! x1 and x2 should be ordered
    if (x(1) > x(2)) then
      call program_error("x not ordered")
    end if

    ! use 1-x2 instead of x2 for better accuracy
    x(2) = 1 - x(2)

    it  = 0
    err = huge(1.0)
    j   = 0.0
    do while ((it < 2) .or. ((err > abs(j) * RTOL) .and. (any(err_dx > abs(x) * RTOL))))
      it = it + 1

      if (it > 20) then
        print "(A,ES25.16E3,A,ES25.16E3,A)", "eta = [", eta(1), ", ", eta(2), "]"
        print "(A,ES25.16E3)", "dpot = ", dpot
        call program_error("No convergence")
      end if

      ! Newton update
      call residual(x, f, dfdx, dfdeta, dfddpot)
      lu = dfdx
      call getrf(lu, ipiv)
      dx = - f
      call getrs(lu, ipiv, dx)
      call gerfs(dfdx, lu, ipiv, -f, dx)
      err_dx = abs(dx)

      ! evaluate current, use error in x to estimate error in j
      call eval_j(x, j, djdx, djddpot)
      err    = sqrt((djdx(1) * dx(1))**2 + (djdx(2) * dx(2))**2)

      ! print "(A,I0,A,ES25.16E3,A,ES25.16E3)", "  it = ", it, "; j = ", j, "; err = ", err
      ! print "(A,2ES25.16E3)", "    x = ", x
      ! print "(A,2ES25.16E3)", "   dx = ", dx

      ! update solution
      xold = x
      x    = x + dx

      ! stay within bounds
      if (x(1) < 0) then
        ! print *, "x1 limited"
        dx = dx * (- xold(1) / dx(1)) * 0.8
        x = xold + dx
      end if
      if (x(2) < 0) then
        ! print *, "x2 limited"
        dx = dx * (- xold(2)) / dx(2) * 0.8
        x = xold + dx
      end if
      if (x(1) > 1 - x(2)) then
        ! print *, "x1,x2 limited"
        dx = dx * (1 - xold(1) - xold(2)) / (dx(1) + dx(2)) * 0.8
        x = xold + dx
      end if
    end do

    ! calculate dxdeta and dxddpot with implicit differentiation
    dxdeta = - dfdeta
    call getrs(lu, ipiv, dxdeta)
    call gerfs(dfdx, lu, ipiv, - dfdeta, dxdeta)
    dxddpot = - dfddpot
    call getrs(lu, ipiv, dxddpot)
    call gerfs(dfdx, lu, ipiv, - dfddpot, dxddpot)

    ! evaluate current
    call eval_j(x, j, djdx, djddpot)
    djdeta(1) =           djdx(1) * dxdeta(1,1) + djdx(2) * dxdeta(2,1)
    djdeta(2) =           djdx(1) * dxdeta(1,2) + djdx(2) * dxdeta(2,2)
    djddpot   = djddpot + djdx(1) * dxddpot(1) + djdx(2) * dxddpot(2)

  contains

    subroutine residual(x, f, dfdx, dfdeta, dfddpot)
      real, intent(in)  :: x(2)
      real, intent(out) :: f(2)
      real, intent(out) :: dfdx(2,2)
      real, intent(out) :: dfdeta(2,2)
      real, intent(out) :: dfddpot(2)

      real          :: dpot1, dpot2, dpot3, j1, j2, j3, dj1deta(2), dj2deta(2), dj3deta(2), dj1ddpot1, dj2ddpot2, dj3ddpot3
      type(hp_real) :: h

      ! get current in left part
      dpot1 = dpot * x(1)
      if (eta(1) < ETA_MICRO) then
        call this%get_REG([eta(1), eta_split(1)], dpot1, j1, dj1deta, dj1ddpot1)
      else
        call this%get_INTERP([eta(1), eta_split(1)], dpot1, j1, dj1deta, dj1ddpot1)
      end if

      ! get current in center part
      h = 1.0 - TwoSum(x(1), x(2))
      dpot2 = dpot * hp_to_real(h)
      call this%get_SGFP([eta_split(1), eta_split(2)], dpot2, j2, dj2deta, dj2ddpot2)

      ! get current in right part
      dpot3 = dpot * x(2)
      if (eta(2) < ETA_MICRO) then
        call this%get_REG([eta_split(2), eta(2)], dpot3, j3, dj3deta, dj3ddpot3)
      else
        call this%get_INTERP([eta_split(2), eta(2)], dpot3, j3, dj3deta, dj3ddpot3)
      end if

      ! residual (respect scale)
      f(1)      = hp_to_real(j1 * h - TwoProduct(j2, x(1)))
      dfdx(1,1) = - j1 - j2 + dj1ddpot1 * dpot2 + dj2ddpot2 * dpot1
      dfdx(1,2) = - j1 + dj2ddpot2 * dpot1
      f(2)      = hp_to_real(j3 * h - TwoProduct(j2, x(2)))
      dfdx(2,1) = - j3 + dj2ddpot2 * dpot3
      dfdx(2,2) = - j2 - j3 + dj2ddpot2 * dpot3 + dj3ddpot3 * dpot2

      ! print "(A,3ES25.16E3)", "    j = ", j1 / x(1), j2 / hp_to_real(h), j3 / x(2)
      ! print "(A,ES25.16E3)",  "   dj = ", abs(j1 / j2 * hp_to_real(h) / x(1) - 1)
      ! print "(A,2ES25.16E3)", "    f = ", f

      dfdeta(1,1) = hp_to_real(dj1deta(1) * h)
      dfdeta(1,2) = 0
      dfdeta(2,1) = 0
      dfdeta(2,2) = hp_to_real(dj3deta(2) * h)

      dfddpot(1) = hp_to_real((dj1ddpot1 - dj2ddpot2) * h * x(1))
      dfddpot(2) = hp_to_real((dj3ddpot3 - dj2ddpot2) * h * x(2))
    end subroutine

    subroutine eval_j(x, j, djdx, djddpot)
      real, intent(in)  :: x(2)
      real, intent(out) :: j
      real, intent(out) :: djdx(2)
      real, intent(out) :: djddpot

      real :: dpot2, j2, dj2deta(2), dj2ddpot2, len

      ! get current in center part
      len = (1.0 - x(1) - x(2))
      dpot2 = dpot * len
      call this%get_SGFP([eta_split(1), eta_split(2)], dpot2, j2, dj2deta, dj2ddpot2)

      j       = j2 / len
      djdx    = (- dj2ddpot2 * dpot + j2 / len) / len
      djddpot = dj2ddpot2
    end subroutine

  end subroutine

  subroutine gauss_table_init(this, min_xi0, max_xi0)
    !! initialize gauss quadrature lookup table
    class(gauss_table), intent(out) :: this
    real,               intent(in)  :: min_xi0
      !! minimal supported xi0
    real,               intent(in)  :: max_xi0
      !! maximal supported xi0

    integer, parameter :: N0 = 1024
    real,    parameter :: RTOL = 1e-15

    integer              :: i, i1, i2, i3, j1, j2, j3, k1, k2, k3, ithread, nthreads, ntot
    integer, allocatable :: nt(:), perm(:)
    real                 :: xi_(NGS), dxi_(NGS), w_(NGS), dw_(NGS), s0_(1), ds0_(1)
    real,    allocatable :: xi0(:), xi(:,:), dxi(:,:), w(:,:), dw(:,:), s0(:), ds0(:)
    type(mpfr_vars)      :: mvars
    type(vector_real)    :: vxi0, vxi, vdxi, vw, vdw, vs0, vds0
    type(vector_int)     :: vi

    ! initial (coarse) grid
    allocate (xi0(N0), xi(NGS,N0), dxi(NGS,N0), w(NGS,N0), dw(NGS,N0), s0(N0), ds0(N0))
    xi0 = logspace(min_xi0, max_xi0, N0)

    !$omp parallel default(none) &
    !$omp private(i,i1,i2,i3,j1,j2,j3,k1,k2,k3,ithread,xi_,dxi_,w_,dw_,s0_,ds0_,mvars,vxi0,vxi,vdxi,vw,vdw,vs0,vds0,vi) &
    !$omp shared(this,nthreads,ntot,nt,xi0,xi,dxi,w,dw,s0,ds0)

    nthreads = omp_get_num_threads()
    ithread  = omp_get_thread_num() + 1

    !$omp single
    allocate (nt(0:nthreads + 1), source = 0)
    !$omp end single

    call mpfr_startup(prec = 256)
    call mvars%init()

    ! generate entries for coarse grid
    !$omp do schedule(dynamic)
    do i = 1, N0
      call this%gen(mvars, xi0(i), xi(:,i), dxi(:,i), w(:,i), dw(:,i), s0(i), ds0(i))
    end do
    !$omp end do

    ! copy coarse grid to thread-local memory
    call vxi0%init(N0,       c = 4 * N0,       x = xi0)
    call vxi%init( N0 * NGS, c = 4 * N0 * NGS, x = reshape(xi,  [N0 * NGS]))
    call vdxi%init(N0 * NGS, c = 4 * N0 * NGS, x = reshape(dxi, [N0 * NGS]))
    call vw%init(  N0 * NGS, c = 4 * N0 * NGS, x = reshape(w,   [N0 * NGS]))
    call vdw%init( N0 * NGS, c = 4 * N0 * NGS, x = reshape(dw,  [N0 * NGS]))
    call vs0%init( N0,       c = 4 * N0,       x = s0)
    call vds0%init(N0,       c = 4 * N0,       x = ds0)

    ! interval stack
    call vi%init(0, c = 16)

    ! refinement
    !$omp do schedule(dynamic)
    do i = 1, N0-1
      ! add interval
      call vi%push(i)
      call vi%push(i+1)

      do while (vi%n > 0)
        ! get interval from stack
        i1   = vi%d(vi%n-1)
        i2   = vi%d(vi%n  )
        vi%n = vi%n - 2

        ! midpoint
        call vxi0%push(exp(0.5*(log(vxi0%d(i1)) + log(vxi0%d(i2)))))
        i3 = vxi0%n

        ! data indices
        j1 = (i1 - 1) * NGS + 1
        j2 = (i2 - 1) * NGS + 1
        j3 = (i3 - 1) * NGS + 1
        k1 = i1 * NGS
        k2 = i2 * NGS
        k3 = i3 * NGS

        ! make room for new point
        call vxi%resize( k3)
        call vdxi%resize(k3)
        call vw%resize(  k3)
        call vdw%resize( k3)
        call vs0%resize( i3)
        call vds0%resize(i3)

        ! generate data for new point and additionally interpolate it from existing data
        call this%gen(mvars, vxi0%d(i3), vxi%d(j3:k3), vdxi%d(j3:k3), vw%d(j3:k3), vdw%d(j3:k3), vs0%d(i3), vds0%d(i3))
        call this%interp(vxi0%d(i1), vxi%d(j1:k1), vdxi%d(j1:k1), vxi0%d(i2), vxi%d(j2:k2), vdxi%d(j2:k2), vxi0%d(i3), xi_, dxi_)
        call this%interp(vxi0%d(i1),  vw%d(j1:k1),  vdw%d(j1:k1), vxi0%d(i2),  vw%d(j2:k2),  vdw%d(j2:k2), vxi0%d(i3),  w_, dw_ )
        call this%interp(vxi0%d(i1), vs0%d(i1:i1), vds0%d(i1:i1), vxi0%d(i2), vs0%d(i2:i2), vds0%d(i2:i2), vxi0%d(i3), s0_, ds0_)

        ! check if refinement is necessary
        if (any(abs(xi_ - vxi%d(j3:k3)) / vxi%d(j3:k3) > RTOL) .or. any(abs(w_ - vw%d(j3:k3)) / vw%d(j3:k3) > RTOL) .or. (abs(s0_(1) - vs0%d(i3)) / vs0%d(i3) > RTOL)) then
          ! add 2 new intervals to stack
          call vi%push(i1)
          call vi%push(i3)
          call vi%push(i3)
          call vi%push(i2)
        end if
      end do
    end do
    !$omp end do

    nt(ithread+1) = vxi0%n - N0
    !$omp barrier

    !$omp single
    ! count elements
    nt(0) = 1
    nt(1) = N0
    do i = 1, nthreads + 1
      nt(i) = nt(i) + nt(i-1)
    end do
    ntot = nt(nthreads + 1) - 1

    ! allocate global memory + fill in coarse grid
    allocate (this%xi0(ntot), this%xi(NGS,ntot), this%dxi(NGS,ntot), this%w(NGS,ntot), this%dw(NGS,ntot), this%s0(ntot), this%ds0(ntot))
    this%xi0(  1:N0) = xi0
    this%xi( :,1:N0) = xi
    this%dxi(:,1:N0) = dxi
    this%w(  :,1:N0) = w
    this%dw( :,1:N0) = dw
    this%s0(   1:N0) = s0
    this%ds0(  1:N0) = ds0
    !$omp end single

    ! copy local values to global memory
    j1 = N0 * NGS + 1
    k1 = vxi0%n * NGS
    this%xi0(  nt(ithread):nt(ithread+1)-1) = vxi0%d(N0+1:vxi0%n)
    this%xi( :,nt(ithread):nt(ithread+1)-1) = reshape(vxi%d( j1:k1), [NGS, vxi0%n - N0])
    this%dxi(:,nt(ithread):nt(ithread+1)-1) = reshape(vdxi%d(j1:k1), [NGS, vxi0%n - N0])
    this%w(  :,nt(ithread):nt(ithread+1)-1) = reshape(vw%d(  j1:k1), [NGS, vxi0%n - N0])
    this%dw( :,nt(ithread):nt(ithread+1)-1) = reshape(vdw%d( j1:k1), [NGS, vxi0%n - N0])
    this%s0(   nt(ithread):nt(ithread+1)-1) = vs0%d( N0+1:vxi0%n)
    this%ds0(  nt(ithread):nt(ithread+1)-1) = vds0%d(N0+1:vxi0%n)

    ! cleanup thread-local memory
    call vxi0%destruct()
    call vxi%destruct()
    call vdxi%destruct()
    call vw%destruct()
    call vdw%destruct()
    call vs0%destruct()
    call vds0%destruct()
    call mvars%destruct()
    call mpfr_cleanup()

    !$omp end parallel

    ! sort
    allocate (perm(size(this%xi0)))
    call qsort(this%xi0, perm = perm)
    this%xi  = this%xi( :,perm)
    this%dxi = this%dxi(:,perm)
    this%w   = this%w(  :,perm)
    this%dw  = this%dw( :,perm)
    this%s0  = this%s0(   perm)
    this%ds0 = this%ds0(  perm)
  end subroutine

  subroutine gauss_table_gen(this, mvars, xi0, xi, dxi, w, dw, s0, ds0)
    !! generate entry
    class(gauss_table), intent(in)    :: this
    type(mpfr_vars),    intent(inout) :: mvars
      !! MPFR temporary values
    real,               intent(in)    :: xi0
      !! weight function parameter
    real,               intent(out)   :: xi(:)
      !! output gauss node positions
    real,               intent(out)   :: dxi(:)
      !! output derivatives of xi wrt xi0
    real,               intent(out)   :: w(:)
      !! output gauss weigmts
    real,               intent(out)   :: dw(:)
      !! output derivatives of w wrt xi0
    real,               intent(out)   :: s0
      !! output 0-th moment (scaling for w)
    real,               intent(out)   :: ds0
      !! output derivative of s0 wrt xi0

    m4_ignore(this)

    ! calculate moments using multiprecision
    call moments()

    ! get recurrence relations of orthogonal polynomials
    call recurrence()

    ! solve tridiagonal eigenvalue system
    call eigenvalues()

  contains

    subroutine moments()
      integer, parameter :: NP = 64

      integer :: m, n

      ! s(0) = log1p(1/xi0)
      ! ds(0) = - 1 / (xi0 * (xi0 + 1))
      call mvars%xi%set(xi0)
      call div(mvars%ix, 1, mvars%xi)
      call log1p_mpfr(mvars%s(0), mvars%ix)
      call add(mvars%t, mvars%xi, 1)
      call mul(mvars%t, mvars%t, xi0)
      call div(mvars%ds(0), -1, mvars%t)

      if (xi0 <= 2.0) then
        ! xi = (-xi0)**n
        ! ix = (-xi0)**(-n)
        ! ss, dss: sum up to this point, derivative
        ! s, ds: (-xi0)**n * ss, derivative
        call mvars%ss%set(mvars%s(0))
        call mvars%dss%set(mvars%ds(0))
        call neg(mvars%xi, mvars%xi)
        call neg(mvars%ix, mvars%ix)

        ! n = 1
        call add(mvars%ss, mvars%ss, mvars%ix)
        call div(mvars%ix, mvars%ix, -xi0)
        call add(mvars%dss, mvars%dss, mvars%ix)
        call mul(mvars%s(1), mvars%xi, mvars%ss)
        call mul(mvars%t, mvars%xi, mvars%dss)
        call sub(mvars%ds(1), mvars%t, mvars%ss)

        ! rest
        do n = 2, MAXN
          call div(mvars%t, mvars%ix, n)
          call add(mvars%ss, mvars%ss, mvars%t)
          call div(mvars%ix, mvars%ix, -xi0)
          call add(mvars%dss, mvars%dss, mvars%ix)
          call mul(mvars%ds(n), mvars%xi, mvars%ss)
          call mul(mvars%ds(n), mvars%ds(n), -n)
          call mul(mvars%xi, mvars%xi, -xi0)
          call mul(mvars%t, mvars%xi, mvars%dss)
          call add(mvars%ds(n), mvars%ds(n), mvars%t)
          call mul(mvars%s(n), mvars%xi, mvars%ss)
        end do
      else
        ! "Convergence Acceleration of Alternating Series", Cohen et al.
        call mvars%d%set(8)
        call sqrt_mpfr(mvars%d, mvars%d)
        call add(mvars%d, mvars%d, 3)
        call pow(mvars%d, mvars%d, NP)
        call div(mvars%t, 1, mvars%d)
        call add(mvars%d, mvars%d, mvars%t)
        call mul(mvars%d, mvars%d, 0.5)
        call mvars%b%set(-1)
        call neg(mvars%c, mvars%d)
        do n = 1, MAXN
          call mvars%s(n)%set(0)
          call mvars%ds(n)%set(0)
        end do
        do m = 0, NP - 1
          call sub(mvars%c, mvars%b, mvars%c)
          do n = 1, MAXN
            call div(mvars%t, mvars%c, n + m + 1)               ! t = c / (n + m + 1)
            call fma(mvars%s(n), mvars%t, mvars%ix, mvars%s(n)) ! s(n) = s(n) + c * ix / (n + m + 1)
          end do
          call div(mvars%ix, mvars%ix, xi0)
          do n = 1, MAXN
            call mul(mvars%t, mvars%c, - (m + 1))                 ! t = - c * (m + 1)
            call div(mvars%t, mvars%t, n + m + 1)                 ! t = - c * (m + 1) / (n + m + 1)
            call fma(mvars%ds(n), mvars%t, mvars%ix, mvars%ds(n)) ! ds(n) = ds(n) - c * ix * (m + 1) / (n + m + 1)
          end do
          call mul(mvars%b, mvars%b, (m + NP) * (m - NP))
          call div(mvars%b, mvars%b, (m + 0.5) * (m + 1))
        end do
        do n = 1, MAXN
          call div( mvars%s(n),  mvars%s(n), mvars%d)
          call div(mvars%ds(n), mvars%ds(n), mvars%d)
        end do
      end if
    end subroutine

    subroutine recurrence()
      integer :: i, j, k

      do i = 1, NGS
        do j = i, NGS + 1
          call mvars%sg( i,j)%set(0)
          call mvars%dsg(i,j)%set(0)
          do k = 1, i - 1
            call mul(mvars%t, mvars%sg(k,i), mvars%sg(k,j))                            ! t = sg(k,i) * sg(k,j)
            call div(mvars%t, mvars%t, mvars%sg(k,k))                                  ! t = sg(k,i) * sg(k,j) / sg(k,k)
            call add(mvars%sg(i,j), mvars%sg(i,j), mvars%t)                            ! sg(i,j) = sg(i,j) + sg(k,i) * sg(k,j) / sg(k,k)
            call fmms(mvars%t, mvars%sg(k,i), mvars%dsg(k,j), mvars%t, mvars%dsg(k,k)) ! t = sg(k,i) * dsg(k,j) - sg(k,i) * sg(k,j) / sg(k,k) * dsg(k,k)
            call fma(mvars%t, mvars%dsg(k,i), mvars%sg(k,j), mvars%t)                  ! t = dsg(k,i) * sg(k,j) + sg(k,i) * dsg(k,j) - sg(k,i) * sg(k,j) / sg(k,k) * dsg(k,k)
            call div(mvars%t, mvars%t, mvars%sg(k,k))                                  ! t = (dsg(k,i) * sg(k,j) + sg(k,i) * dsg(k,j) - sg(k,i) * sg(k,j) / sg(k,k) * dsg(k,k)) / sg(k,k)
            call add(mvars%dsg(i,j), mvars%dsg(i,j), mvars%t)                          ! dsg(i,j) = dsg(i,j) + (dsg(k,i) * sg(k,j) + sg(k,i) * dsg(k,j) - sg(k,i) * sg(k,j) / sg(k,k) * dsg(k,k)) / sg(k,k)
          end do

          ! sg(i,j) = s(i+j-2) - sg(i,j)
          call sub(mvars%sg(i,j), mvars%s(i+j-2), mvars%sg(i,j))

          ! dsg(i,j) = ds(i+j-2) - dsg(i,j)
          call sub(mvars%dsg(i,j), mvars%ds(i+j-2), mvars%dsg(i,j))
        end do
      end do

      call mvars%t%set(0)
      call mvars%u%set(0)
      do i = 1, NGS
        ! c_i = sg(i,i+1) / sg(i,i)
        call div(mvars%c, mvars%sg(i,i+1), mvars%sg(i,i))

        ! dc_i = (dsg(i,i+1) - c * dsg(i,i)) / sg(i,i)
        call mul(mvars%dc, mvars%c, mvars%dsg(i,i))
        call sub(mvars%dc, mvars%dsg(i,i+1), mvars%dc)
        call div(mvars%dc, mvars%dc, mvars%sg(i,i))

        !  a(i) =  c_i -  c_{i-1} =  c - t
        ! da(i) = dc_i - dc_{i-1} = dc - u
        call add(mvars%t, mvars%t,  mvars%c)
        call add(mvars%u, mvars%u, mvars%dc)
        call mvars%aa(i)%set(mvars%t)
        call mvars%daa(i)%set(mvars%u)
        call neg(mvars%t,  mvars%c)
        call neg(mvars%u, mvars%dc)
      end do

      do i = 1, NGS-1
        !  c_i = sg(i+1,i+1) / sg(i,i)
        call div(mvars%t, 1.0, mvars%sg(i,i))                           ! t = 1 / sg(i,i)
        call mul(mvars%c, mvars%sg(i+1,i+1), mvars%t)                   ! c = sg(i+1,i+1) / sg(i,i)

        ! dc_i = (dsg(i+1,i+1) - c_i * dsg(i,i)) / sg(i,i)
        call fms(mvars%dc, mvars%c, mvars%dsg(i,i), mvars%dsg(i+1,i+1)) ! dc = sg(i+1,i+1) / sg(i,i) * dsg(i,i) - dsg(i+1,i+1)
        call mul(mvars%dc, mvars%dc, mvars%t)                           ! dc = (sg(i+1,i+1) / sg(i,i) * dsg(i,i) - dsg(i+1,i+1)) / sg(i,i)
        call neg(mvars%dc, mvars%dc)                                    ! dc = (dsg(i+1,i+1) - sg(i+1,i+1) / sg(i,i) * dsg(i,i)) / sg(i,i)

        !  b_i = sqrt(c_i)
        ! db_i = 1 / (2 * b_i) * dc_i
        call sqrt_mpfr(mvars%t, mvars%c)
        call add(mvars%u, mvars%t, mvars%t)
        call div(mvars%u, mvars%dc, mvars%u)
        call mvars%bb(i)%set(mvars%t)
        call mvars%dbb(i)%set(mvars%u)
      end do
    end subroutine

    subroutine eigenvalues()
      real, parameter :: TOL = 1e-16

      integer :: i, j
      real    :: da(NGS), b(NGS-1), db(NGS-1), v(NGS,NGS), dAv(NGS,NGS), vdAv(NGS,NGS)

      ! approximate matrix
      do i = 1, NGS-1
        xi(i) = mvars%aa( i)%to_real()
        da(i) = mvars%daa(i)%to_real()
        b( i) = mvars%bb( i)%to_real()
        db(i) = mvars%dbb(i)%to_real()
      end do
      xi(NGS) = mvars%aa( NGS)%to_real()
      da(NGS) = mvars%daa(NGS)%to_real()

      ! solve symmetric tridiagonal eigenvalue problem with LAPACK (xi and b are overwritten)
      call stev(xi, b, v)

      ! refine using single Rayleigh quotient iteration step (multiprecision)
      do j = 1, NGS
        ! load estimate
        call mvars%xi%set(xi(j))
        do i = 1, NGS
          call mvars%v(i)%set(v(i,j))
        end do

        ! forward substitution
        call sub(mvars%t, mvars%aa(1), mvars%xi)    ! t = a(1) - xi
        call div(mvars%t, 1.0, mvars%t)             ! t = 1 / (a(1) - xi)
        call mul(mvars%cc(1), mvars%bb(1), mvars%t) ! c(1) = b(1) / (a(1) - xi)
        call mul(mvars%v( 1), mvars%v( 1), mvars%t) ! v(1) = v(1) / (a(1) - xi)
        do i = 2, NGS-1
          call sub(mvars%t, mvars%aa(i), mvars%xi)                   ! t = a(i) - xi
          call fms(mvars%t, mvars%bb(i-1), mvars%cc(i-1), mvars%t)   ! t = b(i-1) * c(i-1) - (a(i) - xi)
          call div(mvars%t, 1.0, mvars%t)                            ! t = 1 / (b(i-1) * c(i-1) - (a(i) - xi))
          call mul(mvars%cc(i), mvars%bb(i), mvars%t)                ! c(i) = b(i) / (b(i-1) * c(i-1) - (a(i) - xi))
          call neg(mvars%cc(i), mvars%cc(i))                         ! c(i) = b(i) / (a(i) - xi - b(i-1) * c(i-1))
          call fms(mvars%u, mvars%bb(i-1), mvars%v(i-1), mvars%v(i)) ! u = b(i-1) * v(i-1) - v(i)
          call mul(mvars%v(i), mvars%u, mvars%t)                     ! v(i) = (v(i) - b(i-1) * v(i-1)) / (a(i) - xi - b(i-1) * c(i-1))
        end do
        call sub(mvars%t, mvars%aa(NGS), mvars%xi)                       ! t = a(N) - xi
        call fms(mvars%t, mvars%bb(NGS-1), mvars%cc(NGS-1), mvars%t)     ! t = b(N-1) * c(N-1) - (a(N) - xi)
        call fms(mvars%u, mvars%bb(NGS-1), mvars%v(NGS-1), mvars%v(NGS)) ! u = b(N-1) * v(N-1) - v(N)
        call div(mvars%v(NGS), mvars%u, mvars%t)                         ! v(N) = (v(N) - b(N-1) * v(N-1)) / (a(N) - xi - b(N-1) * c(N-1))

        ! back substitution
        do i = NGS-1, 1, -1
          call fms(mvars%v(i), mvars%cc(i), mvars%v(i+1), mvars%v(i)) ! v(i) = c(i) * v(i+1) - v(i)
          call neg(mvars%v(i), mvars%v(i))                            ! v(i) = v(i) - c(i) * v(i+1)
        end do

        ! normalization
        call mvars%t%set(0.0)
        do i = 1, NGS
          call sqr(mvars%u, mvars%v(i))
          call add(mvars%t, mvars%t, mvars%u)
        end do
        call sqrt_mpfr(mvars%t, mvars%t)
        do i = 1, NGS
          call div(mvars%v(i), mvars%v(i), mvars%t)
        end do

        ! improve eigenvalue estimate
        call fmma(mvars%t, mvars%aa(1), mvars%v(1), mvars%bb(1), mvars%v(2)) ! t = a(1) * v(1) + b(1) * v(2)
        call mul(mvars%xi, mvars%v(1), mvars%t)                              ! xi = v(1) * (a(1) * v(1) + b(1) * v(2))
        do i = 2, NGS-1
          call fmma(mvars%t, mvars%bb(i-1), mvars%v(i-1), mvars%aa(i), mvars%v(i)) ! t = b(i-1) * v(i-1) + a(i) * v(i)
          call fma(mvars%t, mvars%bb(i), mvars%v(i+1), mvars%t)                    ! t = b(i-1) * v(i-1) + a(i) * v(i) + b(i) * v(i+1)
          call fma(mvars%xi, mvars%v(i), mvars%t, mvars%xi)                        ! xi = xi + v(i) * (b(i-1) * v(i-1) + a(i) * v(i) + b(i) * v(i+1))
        end do
        call fmma(mvars%t, mvars%bb(NGS-1), mvars%v(NGS-1), mvars%aa(NGS), mvars%v(NGS)) ! t = b(N-1) * v(N-1) + a(N) * v(N)
        call fma(mvars%xi, mvars%v(NGS), mvars%t, mvars%xi)                              ! xi = xi + v(N) * (b(N-1) * v(N-1) + a(N) * v(N))

        ! convert to real
        xi(j) = mvars%xi%to_real()
        do i = 1, NGS
          v(i,j) = mvars%v(i)%to_real()
        end do

        ! ! extract weight (w = v(1)**2 * s(0))
        ! call sqr(mvars%t, mvars%v(1))
        ! call mul(mvars%t, mvars%t, mvars%s(0))
        ! w(j) = mvars%t%to_real()

        ! extract weight w = v(1)**2
        call sqr(mvars%t, mvars%v(1))
        w(j) = mvars%t%to_real()
      end do

      ! dA * v
      do i = 1, NGS
        dAv( :,     i) =                  da * v( :,     i)
        dAv(1:NGS-1,i) = dAv(1:NGS-1,i) + db * v(2:NGS,  i)
        dAv(2:NGS,  i) = dAv(2:NGS,  i) + db * v(1:NGS-1,i)
      end do

      ! v' * dA * v
      call gemm(v, dAv, vdAv, transA = 'T')

      ! dxi and dw
      do i = 1, NGS
        dxi(i) = vdAv(i,i)
        dw(i) = 0
        do j = 1, NGS
          if (j == i) cycle
          dw(i) = dw(i) + vdAv(j,i) / (xi(i) - xi(j)) * v(1,j)
        end do
      end do
      ! dw = v(1,:) * (2 * dw * mvars%s(0)%to_real() + v(1,:) * mvars%ds(0)%to_real())
      dw = 2 * v(1,:) * dw

      ! s0 and ds0
      s0  = mvars%s( 0)%to_real()
      ds0 = mvars%ds(0)%to_real()
    end subroutine

  end subroutine

  subroutine gauss_table_get(this, xi0, xi, dxi, w, dw, s0, ds0)
    class(gauss_table), intent(in)    :: this
    real,               intent(in)    :: xi0
      !! weight function parameter
    real,               intent(out)   :: xi(:)
      !! output gauss node positions
    real,               intent(out)   :: dxi(:)
      !! output derivatives of xi wrt xi0
    real,               intent(out)   :: w(:)
      !! output gauss weigmts
    real,               intent(out)   :: dw(:)
      !! output derivatives of w wrt xi0
    real,               intent(out)   :: s0
      !! output 0-th moment (scaling for w)
    real,               intent(out)   :: ds0
      !! output derivative of s0 wrt xi0

    integer :: k
    real    :: s0_(1), ds0_(1)

    if ((xi0 < this%xi0(1)) .or. (xi0 > this%xi0(size(this%xi0)))) then
      print "(A,ES25.16E3)", "xi0 = ", xi0
      call program_error("xi0 is out of range")
    end if

    k = bin_search(this%xi0, xi0, mode = BS_LESS)

    call this%interp(this%xi0(k), this%xi(:,k), this%dxi(:,k), this%xi0(k+1), this%xi(  :,k+1), this%dxi(  :,k+1), xi0,  xi,  dxi )
    call this%interp(this%xi0(k), this%w( :,k), this%dw( :,k), this%xi0(k+1), this%w(   :,k+1), this%dw(   :,k+1), xi0,   w,  dw  )
    call this%interp(this%xi0(k), this%s0(k:k), this%ds0(k:k), this%xi0(k+1), this%s0(k+1:k+1), this%ds0(k+1:k+1), xi0, s0_,  ds0_)
    s0  = s0_(1)
    ds0 = ds0_(1)
  end subroutine

  subroutine gauss_table_interp(this, x1, f1, df1, x2, f2, df2, x, f, df)
    class(gauss_table), intent(in)  :: this
    real,               intent(in)  :: x1
    real,               intent(in)  :: f1(:)
    real,               intent(in)  :: df1(:)
    real,               intent(in)  :: x2
    real,               intent(in)  :: f2(:)
    real,               intent(in)  :: df2(:)
    real,               intent(in)  :: x
    real,               intent(out) :: f(:)
    real,               intent(out) :: df(:)

    real :: e, e1, e2, de, t
    real :: h00, h10, h01, h11, g00, g10, g01, g11

    m4_ignore(this)

    e  = log(x)
    e1 = log(x1)
    e2 = log(x2)
    de = e2 - e1
    t  = (e - e1) / (e2 - e1)

    h00 = (1 + 2 * t) * (1 - t)**2
    h10 = t * (1 - t)**2
    h01 = t**2 * (3 - 2*t)
    h11 = t**2 * (t - 1)

    g00 = 6 * t * (t - 1)
    g10 = (3 * t - 1) * (t - 1)
    g01 = - 6 * t * (t - 1)
    g11 = t * (3 * t - 2)

    f  = h00 * f1 + h10 * de * x1 * df1 + h01 * f2 + h11 * de * x2 * df2
    df = g00 * f1 + g10 * de * x1 * df1 + g01 * f2 + g11 * de * x2 * df2
    df = df / (de * x)
  end subroutine

  subroutine gauss_table_output(this, fname)
    !! output gauss lookup table to file
    class(gauss_table), intent(in) :: this
    character(*),       intent(in) :: fname

    character(32) :: fmt
    integer       :: i, funit

    write (fmt, "(A,I0,A)") "(", 1 + 4*NGS + 2, "ES25.16E3)"

    open (newunit = funit, file = fname, status = "replace", action = "write")
    do i = 1, size(this%xi0)
      write (funit, fmt) this%xi0(i), this%xi(:,i), this%dxi(:,i), this%w(:,i), this%dw(:,i), this%s0(i), this%ds0(i)
    end do
    close (funit)
  end subroutine

  subroutine mpfr_vars_init(this)
    !! allocate memory for MPFR variables
    class(mpfr_vars), intent(out) :: this

    integer :: i, j, n

    do n = 0, MAXN
      call this%s( n)%init()
      call this%ds(n)%init()
    end do
    do i = 1, NGS
      do j = 1, NGS + 1
        call this%sg( i,j)%init()
        call this%dsg(i,j)%init()
      end do
    end do
    call this%xi%init()
    call this%ix%init()
    call this%ss%init()
    call this%dss%init()
    call this%b%init()
    call this%c%init()
    call this%dc%init()
    call this%d%init()
    call this%t%init()
    call this%u%init()
    do i = 1, NGS
      call this%aa(i)%init()
      call this%v(i)%init()
      call this%daa(i)%init()
    end do
    do i = 1, NGS-1
      call this%bb(i)%init()
      call this%cc(i)%init()
      call this%dbb(i)%init()
    end do
  end subroutine

  subroutine mpfr_vars_destruct(this)
    class(mpfr_vars), intent(inout) :: this

    integer :: i, j, n

    do n = 0, MAXN
      call this%s( n)%destruct()
      call this%ds(n)%destruct()
    end do
    do i = 1, NGS
      do j = 1, NGS + 1
        call this%sg( i,j)%destruct()
        call this%dsg(i,j)%destruct()
      end do
    end do
    call this%xi%destruct()
    call this%ix%destruct()
    call this%ss%destruct()
    call this%dss%destruct()
    call this%b%destruct()
    call this%c%destruct()
    call this%dc%destruct()
    call this%d%destruct()
    call this%t%destruct()
    call this%u%destruct()
    do i = 1, NGS
      call this%aa(i)%destruct()
      call this%v(i)%destruct()
      call this%daa(i)%destruct()
    end do
    do i = 1, NGS-1
      call this%bb(i)%destruct()
      call this%cc(i)%destruct()
      call this%dbb(i)%destruct()
    end do
  end subroutine

end module
