subroutine qsort_{}T{}(x, perm, init_perm)
  m4_type(T),        intent(inout) :: x(:)
    !! array which should be sorted
  integer, optional, intent(inout) :: perm(:)
    !! permutation array
  logical, optional, intent(in)    :: init_perm
    !! initialize permutation array by consecutive numbers (default: true)

  integer :: i
  logical :: init_perm_

  ! init permutation array
  if (present(perm)) then
    m4_assert(size(x) == size(perm))
    init_perm_ = .true.
    if (present(init_perm)) init_perm_ = init_perm
    if (init_perm_) perm = [(i, i=1,size(x))]
  end if

  ! start recursive quick sort
  call qsort_rec(1, size(x))

contains

  recursive subroutine qsort_rec(i0, i1)
    integer, intent(in) :: i0
    integer, intent(in) :: i1

    integer    :: mid, l, r
    m4_type(T) :: pivot

    ! 0 or 1 element
    if (i0 >= i1) return

    ! quicksort with optimal sorting networks and insertion sort as base cases
    select case (i1 - i0)
      case (1) ! 2 elements
        call swap_cmp(i0,   i1  )
      case (2) ! 3 elements
        call swap_cmp(i0+1, i1  )
        call swap_cmp(i0  , i1  )
        call swap_cmp(i0  , i0+1)
      case (3) ! 4 elements
        call swap_cmp(i0  , i0+1)
        call swap_cmp(i1-1, i1  )
        call swap_cmp(i0  , i1-1)
        call swap_cmp(i0+1, i1  )
        call swap_cmp(i0+1, i1-1)
      case (4) ! 5 elements
        call swap_cmp(i0  , i0+1)
        call swap_cmp(i1-1, i1  )
        call swap_cmp(i0+2, i1  )
        call swap_cmp(i0+2, i1-1)
        call swap_cmp(i0  , i1-1)
        call swap_cmp(i0  , i0+2)
        call swap_cmp(i0+1, i1  )
        call swap_cmp(i0+1, i1-1)
        call swap_cmp(i0+1, i0+2)
      case (5:15) ! 6 - 16 elements : insertion sort
        do r = i0 + 1, i1
          pivot = x(r)
          if (present(perm)) mid = perm(r)

          do l = r - 1, i0, -1
            if (x(l) <= pivot) exit
            x(l+1) = x(l)
            if (present(perm)) perm(l+1) = perm(l)
          end do
          x(l+1) = pivot
          if (present(perm)) perm(l+1) = mid
        end do
      case default ! > 16 elements : quick sort
        ! set pivot to median of three (sort them in the process)
        mid = (i0 + i1) / 2
        call swap_cmp(mid, i1 )
        call swap_cmp(i0,  i1 )
        call swap_cmp(i0,  mid)
        pivot = x(mid)

        ! sort elements left or right according to pivot
        l = i0
        r = i1
        do while (.true.)
          l = l + 1
          do while (x(l) < pivot)
            l = l + 1
          end do
          r = r - 1
          do while (x(r) > pivot)
            r = r - 1
          end do
          if (l >= r) exit
          call swap(l, r)
        end do

        ! recursion on two sub-arrays
        call qsort_rec(i0,  r )
        call qsort_rec(r+1, i1)
    end select
  end subroutine

  subroutine swap_cmp(i0, i1)
    integer, intent(in) :: i0
    integer, intent(in) :: i1

    if (x(i1) < x(i0)) call swap(i0, i1)
  end subroutine

  subroutine swap(i0, i1)
    integer, intent(in) :: i0
    integer, intent(in) :: i1

    integer    :: itmp
    m4_type(T) :: tmp

    tmp   = x(i0)
    x(i0) = x(i1)
    x(i1) = tmp
    if (present(perm)) then
      itmp     = perm(i0)
      perm(i0) = perm(i1)
      perm(i1) = itmp
    end if
  end subroutine

end subroutine

m4_undefine({T})
