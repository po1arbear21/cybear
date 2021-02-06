#include "../util/macro.f90.inc"

module simple_equations_m
  use equation_m
  implicit none

  type, extends(equation) :: dummy_equation
    !! dummy equation which provides values but does not change anything (e.g. stat values for dynamic sim)
  contains
    procedure :: init => dummy_equation_init
    procedure :: eval => dummy_equation_eval
  end type

  type, extends(equation) :: selector_equation
    !! equation that selects variables from one or multiple var selectors (e.g. select [n] from [n; j])

    type(dirichlet_stencil) :: st
  contains
    procedure :: init => selector_equation_init
    procedure :: eval => selector_equation_eval
  end type

  type dummy_equation_ptr
    type(dummy_equation), pointer :: p => null()
  end type

  type selector_equation_ptr
    type(selector_equation), pointer :: p => null()
  end type

#define T dummy_equation_ptr
#define TT type(dummy_equation_ptr)
#include "../util/vector_def.f90.inc"

#define T selector_equation_ptr
#define TT type(selector_equation_ptr)
#include "../util/vector_def.f90.inc"

contains

#define T dummy_equation_ptr
#define TT type(dummy_equation_ptr)
#include "../util/vector_imp.f90.inc"

#define T selector_equation_ptr
#define TT type(selector_equation_ptr)
#include "../util/vector_imp.f90.inc"

  subroutine dummy_equation_init(this, v)
    !! initialize dummy equation
    class(dummy_equation),    intent(out) :: this
    class(vselector), target, intent(in)  :: v
      !! dummy variable

    ! init base
    call this%equation_init("Provide"//v%name)

    ! add provided var
    call this%add_prov(v)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine dummy_equation_eval(this)
    !! evaluate dummy equation
    class(dummy_equation), intent(inout) :: this

    ! do nothing
    IGNORE(this)
  end subroutine

  subroutine selector_equation_init(this, v1, v2, status)
    !! initialize selector equation
    class(selector_equation), target, intent(out) :: this
    type(vselector),          target, intent(in)  :: v1
      !! result var selector (select v1 ...)
    type(vselector_ptr),              intent(in)  :: v2(:)
      !! dependency var selectors (... from v2(:))
    logical,                          intent(out) :: status
      !! return success/fail: success=true, fail=false

    integer           :: i, j, ival1, ival2, itab1, itab2(v1%ntab)
    type(stencil_ptr) :: st(v1%ntab)

    ! init base
    call this%equation_init("Select"//v1%name)

    ! provide v1
    call this%add_prov(v1)

    ! stencil
    call this%st%init() ! dirichlet stencil
    do itab1 = 1, v1%ntab
      st(itab1)%p => this%st
    end do

    ! loop over v1 variables
    do ival1 = 1, v1%nval
      ! loop over source var selectors
      do i = 1, size(v2)
        block
          logical :: valmsk(v1%nval,v2(i)%p%nval)

          ! find variable in v2(i)%p%v
          ival2 = -1
          do j = 1, v2(i)%p%nval
            if (associated(v1%v(ival1)%p, v2(i)%p%v(j)%p)) then
              ival2 = j
              exit
            end if
          end do
          if (ival2 <= 0) cycle

          ! check tables
          itab2 = -1
          do itab1 = 1, v1%ntab
            do j = 1, v2(i)%p%ntab
              if (associated(v1%tab(itab1)%p, v2(i)%p%tab(j)%p)) then
                itab2(itab1) = j
                exit
              end if
            end do
          end do
          if (any(itab2 <= 0)) cycle

          ! add dependency
          call this%add_dep(v2(i)%p)

          ! init jacobian
          valmsk = .false.
          valmsk(ival1,ival2) = .true.
          call this%init_jaco(1, this%vdep%n, st, const = .true., valmsk = valmsk)

          ! set jacobian entries
          do itab1 = 1, v1%ntab
            do j = 1, v1%tab(itab1)%p%n
              ! set derivative to 1
              this%jaco(1,this%vdep%n)%p%d(itab1)%d(ival1,ival2,1,j) = 1.0
            end do
          end do
        end block
      end do

      if (i > size(v2)) then
        status = .false.
        call this%destruct() ! clean up
        return
      end if
    end do

    ! finish initialization
    status = .true.
    call this%init_final()
  end subroutine

  subroutine selector_equation_eval(this)
    !! evaluate selector equation
    class(selector_equation), intent(inout) :: this

    ! do nothing
    IGNORE(this)
  end subroutine

end module