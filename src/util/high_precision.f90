module high_precision_m
  !! algorithms copied from "Accurate sum and dot prodcut" by ogita, rump and oishi

  use error_m

  implicit none

  private
  public :: hp_dot
  public :: hp_sum
  public :: hp_to_real
  public :: real_to_hp
  public :: operator(+), operator(-)

  public :: hp_real

  type hp_real
    !! represents high precision value x + c

    real :: x
      !! principal value
    real :: y
      !! correction value
  end type

  interface operator (+)
    module procedure :: hp_real_add_hh
    module procedure :: hp_real_add_hr
    module procedure :: hp_real_add_rh
  end interface

  interface operator (-)
    module procedure :: hp_real_neg
    module procedure :: hp_real_sub_hh
    module procedure :: hp_real_sub_hr
    module procedure :: hp_real_sub_rh
  end interface

  interface hp_sum
    module procedure :: SumKvert
  end interface

  interface hp_dot
    module procedure :: DotK
  end interface

contains

  elemental function TwoSum(a, b) result(h)
    !! error free transformation of the sum of two floating point numbers

    real, intent(in) :: a
    real, intent(in) :: b
    type(hp_real)    :: h

    real             :: z

    ! principal
    h%x = a + b

    ! correction
    z   = h%x - a
    h%y = (a - (h%x - z)) + (b - z)
  end function

  elemental function Split(a) result(h)
    !! error free splotting of float into two parts

    real, intent(in) :: a
    type(hp_real)    :: h

    real            :: c
    real, parameter :: factor = 134217729.0      ! 2^27+1

    c   = factor * a
    h%x = c-(c-a)
    h%y = a-h%x
  end function

  elemental function TwoProduct(a, b) result(h)
    !! multiply two real numbers with high precision (error free product)

    real, intent(in) :: a
    real, intent(in) :: b
    type(hp_real)    :: h

    type(hp_real) :: h1, h2

    ! principal
    h%x = a * b

    ! correction
    h1 = Split(a)
    h2 = Split(b)
    associate(y => h%y, x=> h%x, a1 => h1%x, a2 => h1%y, b1 => h2%x, b2 => h2%y)
      y = a2*b2 - (((x-a1*b1) - a2*b1) - a1*b2)
    end associate
  end function

  function SumKvert(p, K) result(res)
    !! high precision K-fold summation.

    real              , intent(in)  :: p(:)
    integer, optional , intent(in)  :: K
      !! K-fold precision.
      !! default: 2
    real                            :: res

    integer           :: n, i, j, k_, KK
    real              :: s, alp
    real, allocatable :: q(:)
    type(hp_real)     :: h_tmp

    KK = 2
    if (present(K)) KK = K

    if (KK <  1) then
      call program_error('K should be greater than 0')
    else if (KK == 1) then
      res = sum(p)
      return
    end if

    n   = size(p)
    KK  = min(KK, n)
    allocate(q(KK-1))

    do i = 1, KK-1
      s = p(i)
      do k_ = 1, i-1
        ! [qk, s] <- TowSum(qk, s)
        h_tmp = TwoSum(q(k_), s)
        q(k_) = h_tmp%x
        s     = h_tmp%y
      end do
      q(i) = s
    end do

    s = 0.0
    do i = KK, n
      alp = p(i)
      do k_ = 1, KK-1
        ! [qk, alpha] <- TowSum(qk, alpha)
        h_tmp = TwoSum(q(k_), alp)
        q(k_) = h_tmp%x
        alp   = h_tmp%y
      end do
      s = s + alp
    end do

    do j = 1, KK-2
      alp = q(j)
      do k_ = j+1, KK-1
        ! [qk, alpha] <- TowSum(qk, alpha)
        h_tmp = TwoSum(q(k_), alp)
        q(k_) = h_tmp%x
        alp   = h_tmp%y
      end do
      s = s + alp
    end do

    res = s + q(KK-1)
  end function

  function DotK(x, y, K) result(res)
    !! high precision K-fold dot product.

    real              , intent(in)  :: x(:)
    real              , intent(in)  :: y(:)
    integer, optional , intent(in)  :: K
      !! K-fold precision.
      !! default: 3
    real                            :: res

    integer           :: n, i, K_
    real              :: p, h
    real, allocatable :: r(:)
    type(hp_real)     :: hp_tmp

    K_ = 3
    if (present(K)) K_ = K

    if (K_ <  3) call program_error('K must be greater than 2')

    n = size(x)
    allocate(r(2*n))

    ! [p, r1] = TwoProduct(x1, y1)
    hp_tmp  = TwoProduct(x(1), y(1))
    p       = hp_tmp%x
    r(1)    = hp_tmp%y

    do i = 2, n
      ! [h, ri] = TwoProduct(xi, yi)
      hp_tmp  = TwoProduct(x(i), y(i))
      h       = hp_tmp%x
      r(i)    = hp_tmp%y

      ! [p, r_{n+i-1}] = TwoSum(p, h)
      hp_tmp    = TwoSum(p, h)
      p         = hp_tmp%x
      r(n+i-1)  = hp_tmp%y
    end do
    r(2*n) = p
    res = SumKvert(r, K_-1)
  end function



  elemental function real_to_hp(r) result(h)
    !! convert real to high precision real

    real, intent(in) :: r
    type(hp_real)    :: h

    ! principal is r; correction zero
    h%x = r
    h%y = 0
  end function

  elemental function hp_to_real(h) result(r)
    !! convert high precision real to real

    type(hp_real), intent(in) :: h
    real                      :: r

    ! add principal and correction
    r = h%x + h%y
  end function

  elemental function hp_real_add_hh(h1, h2) result(h3)
    !! add two high precision reals

    type(hp_real), intent(in) :: h1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3    = h1 + h2%x       ! calls around hp_real_add_hr
    h3%y  = h3%y + h2%y
  end function

  elemental function hp_real_add_hr(h1, r2) result(h3)
    !! add real to high precision real

    type(hp_real), intent(in) :: h1
    real,          intent(in) :: r2
    type(hp_real)             :: h3

    h3 = TwoSum(h1%x, r2)

    ! update correction
    h3%y = h3%y + h1%y
  end function

  elemental function hp_real_add_rh(r1, h2) result(h3)
    !! add high precision real to real

    real,          intent(in) :: r1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3


    h3 = h2 + r1          ! calls hp_real_add_hr
  end function

  elemental function hp_real_neg(h1) result(h2)
    !! negate high precision real

    type(hp_real), intent(in) :: h1
    type(hp_real)             :: h2

    ! negate principal and correction
    h2%x = - h1%x
    h2%y = - h1%y
  end function

  elemental function hp_real_sub_hh(h1, h2) result(h3)
    !! subtract two high precision reals

    type(hp_real), intent(in) :: h1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = h1 + (-h2)           ! calls hp_real_neg -> hp_real_add_hh
  end function

  elemental function hp_real_sub_hr(h1, r2) result(h3)
    !! subtract real from high precision real

    type(hp_real), intent(in) :: h1
    real,          intent(in) :: r2
    type(hp_real)             :: h3

    h3 = h1 + (-r2)           ! calls hp_real_add_hr
  end function

  elemental function hp_real_sub_rh(r1, h2) result(h3)
    !! subtract high precision real from real

    real,          intent(in) :: r1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = r1 + (-h2)           ! calls hp_real_neg -> hp_real_add_rh
  end function

end module high_precision_m
