module jacobian_chain_m
  use jacobian_matrix_m
  implicit none

  type, abstract :: jacobian_chain
    !! applies the chain-rule between jacobians (in matrix representation)
    type(jacobian_matrix) :: result
      !! resulting jacobian matrix
  contains
    procedure :: init => jacobian_chain_init
    procedure :: eval => jacobian_chain_eval
  end type

  type, extends(jacobian_chain) :: jacobian_add_chain
    !! addition
    type(jacobian_matrix_ptr), allocatable :: jaco(:)
      !! summands
  end type

  type, extends(jacobian_chain) :: jacobian_mul_chain
    !! multiplication
    type(jacobian_matrix), pointer :: jaco1
      !! first matrix
    type(jacobian_matrix), pointer :: jaco2
      !! second matrix
  end type

  type jacobian_chain_ptr
    !! polymorphic pointer to jacobian chain
    class(jacobian_chain), pointer :: p
  end type

#define T jacobian_chain_ptr
#define TT type(jacobian_chain_ptr)
#include "../util/vector_def.f90.inc"

contains

#define T jacobian_chain_ptr
#define TT type(jacobian_chain_ptr)
#include "../util/vector_imp.f90.inc"

  subroutine jacobian_chain_init(this, jaco)
    !! initialize jacobian chain
    class(jacobian_chain),     intent(out) :: this
    type(jacobian_matrix_ptr), intent(in)  :: jaco(:)
      !! input jacobian matrices

    integer :: i, ival1, ival2, ival3, itab1, itab2, itab3, nval1, nval2, ntab1, ntab2

    nval1 = jaco(1)%p%v1%nval
    nval2 = jaco(2)%p%v2%nval
    ntab1 = jaco(1)%p%v1%ntab
    ntab2 = jaco(2)%p%v2%ntab

    block
      logical :: const(ntab1,ntab2), zero(ntab1,ntab2), valmsk(nval1,nval2,ntab1,ntab2)

      const  = .true.
      zero   = .true.
      valmsk = .false.

      select type (this)
        type is (jacobian_add_chain)
          ! save pointers to input jacobian matrices
          this%jaco = jaco

          ! result is const/zero if all summands are const/zero; values exist if any of summands has values
          do i = 1, size(jaco)
            const  = const  .and. jaco(i)%p%const
            zero   = zero   .and. jaco(i)%p%zero
            valmsk = valmsk .or.  jaco(i)%p%valmsk
          end do

        type is (jacobian_mul_chain)
          ! save pointers to input jacobian matrices
          this%jaco1 => jaco(1)%p
          this%jaco2 => jaco(2)%p

          ! perform pseudo matrix multiplication to figure out const, zero and valmsk flags
          do itab1 = 1, ntab1; do itab3 = 1, jaco(1)%p%v2%ntab
            if (jaco(1)%p%zero(itab1,itab3)) cycle
            do itab2 = 1, ntab2
              if (jaco(2)%p%zero(itab3,itab2)) cycle
              const(itab1,itab2) = const(itab1,itab2) &
                &  .and. jaco(1)%p%const(itab1,itab3) &
                &  .and. jaco(2)%p%const(itab3,itab2)
              zero(itab1,itab2) = .false.
              do ival1 = 1, jaco(1)%p%v1%nval; do ival3 = 1, jaco(1)%p%v2%nval
                do ival2 = 1, jaco(2)%p%v2%nval
                  valmsk(ival1,ival2,itab1,itab2) = valmsk(ival1,ival2,itab1,itab2) .or. &
                    &                    (jaco(1)%p%valmsk(ival1,ival3,itab1,itab3) .and. &
                    &                     jaco(2)%p%valmsk(ival3,ival2,itab3,itab2))
                end do
              end do; end do
            end do
          end do; end do
      end select

      ! initialize result matrix
      call this%result%init(jaco(1)%p%v1, jaco(2)%p%v2, const, zero, valmsk)

      ! set constant parts of result matrix
      call this%eval(nonconst = .false.)
    end block
  end subroutine

  subroutine jacobian_chain_eval(this, const, nonconst)
    !! evaluate jacobian chain
    class(jacobian_chain), intent(inout) :: this
    logical, optional,     intent(in)    :: const
      !! evaluate const blocks (default: true)
    logical, optional,     intent(in)    :: nonconst
      !! evaluate non-const blocks (default: true)

    integer :: itab1, itab2
    logical :: const_, nonconst_, ev(this%result%v1%ntab,this%result%v2%ntab)

    ! optional arguments
    const_ = .true.
    if (present(const)) const_ = const
    nonconst_ = .true.
    if (present(nonconst)) nonconst_ = nonconst

    ! decide which blocks to evaluate
    ev = (.not. this%result%zero) .and. ((this%result%const .and. const_) .or. (.not. this%result%const .and. nonconst_))

    ! evaluate
    select type (this)
      type is (jacobian_add_chain)
        ! add matrices
        do itab1 = 1, this%result%v1%ntab; do itab2 = 1, this%result%v2%ntab
          if (ev(itab1,itab2)) call this%result%add_jaco(this%jaco, itab1 = itab1, itab2 = itab2)
        end do; end do
      type is (jacobian_mul_chain)
        ! multiply matrices
        do itab1 = 1, this%result%v1%ntab; do itab2 = 1, this%result%v2%ntab
          if (ev(itab1,itab2)) call this%result%mul_jaco(this%jaco1, this%jaco2, itab1 = itab1, itab2 = itab2)
        end do; end do
    end select
  end subroutine

end module