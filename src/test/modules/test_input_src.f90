module test_input_src_m

  use input_src_m
  use math_m,      only: linspace, PI
  use test_case_m, only: test_case

  implicit none

  private
  public test_input_src

contains

  subroutine test_input_src()
    type(test_case) :: tc

    call tc%init("input_src")

    ! const_src
    block
      type(const_src) :: src

      ! test init & get
      call src%init([5.0])

      call tc%assert_eq(src%get(10.0), [5.0], 0.0, 0.0, "const: get 1")
      call tc%assert_eq(src%get(-1.0), [5.0], 0.0, 0.0, "const: get 2")

      ! test reset with init
      call src%init([-2.0])
      call tc%assert_eq(src%get(10.0), [-2.0], 0.0, 0.0, "const: init/reset")
    end block

    ! polygon_src
    block
      real              :: t(5), y(1,5)
      type(polygon_src) :: src

      t      = linspace(0.0, 4.0, 5)
      y(1,:) = linspace(0.0, 8.0, 5)

      ! test init & get
      call src%init(t, y)

      call tc%assert_eq(src%get(0.0), [0.0], 0.0, 0.0, "polygon: get 1")
      call tc%assert_eq(src%get(2.0), [4.0], 0.0, 0.0, "polygon: get 2")
      call tc%assert_eq(src%get(3.5), [7.0], 0.0, 0.0, "polygon: get 3")
    end block

    ! sine_src
    block
      real           :: freq, ampl(2), phase(2)
      type(sine_src) :: src

      freq  = 1.0
      ampl  = [1.0, -2.0]
      phase = [0.0, PI]

      ! test init & get
      call src%init(freq, ampl, phase)

      call tc%assert_eq(src%get(0.0),  [ 0.0,  0.0], 1e-14, 3e-16, "sine: get 1")
      call tc%assert_eq(src%get(0.25), [ 1.0,  2.0], 1e-14, 3e-16, "sine: get 2")
      call tc%assert_eq(src%get(0.75), [-1.0, -2.0], 1e-14, 3e-16, "sine: get 3")
    end block

    ! harmonic_src
    block
      real               :: freq, c(1,0:2), s(1,2)
      type(harmonic_src) :: src

      freq  = 1.0
      c     = reshape([1.0,  -2.0, 0.5], [1,3])
      s     = reshape(      [ 2.0, 1.0], [1,2])

      ! test init & get
      call src%init(freq, c, s)

      call tc%assert_eq(src%get(0.0),  [-0.5], 1e-14, 1e-16, "harmonic: get 1")
      call tc%assert_eq(src%get(0.25), [ 2.5], 1e-14, 1e-16, "harmonic: get 2")
      call tc%assert_eq(src%get(0.75), [-1.5], 1e-14, 1e-16, "harmonic: get 3")
    end block

    ! periodic_polygon
    block
      real                       :: freq, tn(3), y(1,3)
      type(periodic_polygon_src) :: src

      ! check reasonable values
      freq = 1.0
      tn   = [0.0, 0.25, 0.75]/freq
      y    = reshape([0.0, 5.0, 10.0], [1,3])

      ! test init & get
      call src%init(freq, tn, y)

      ! testing these values
      call tc%assert_eq(src%get(0.0), [0.0], 1e-14, 1e-16, "periodic_polygon: get 1")
      call tc%assert_eq(src%get(0.1), [2.0], 1e-14, 1e-16, "periodic_polygon: get 2")
      call tc%assert_eq(src%get(0.5), [7.5], 1e-14, 1e-16, "periodic_polygon: get 3")
      call tc%assert_eq(src%get(1.0), [0.0], 1e-14, 1e-16, "periodic_polygon: get 4")
    end block

    call tc%finish()
  end subroutine

end module
