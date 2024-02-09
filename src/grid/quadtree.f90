submodule(triang_grid_m) quadtree_sm

  use math_m,    only: cross_product_2d
  use vector_m,  only: vector_int

  implicit none

  ! for better readability:
  integer, parameter :: LOWER = 1, UPPER = 2
  integer, parameter :: X_IDX = 1, Y_IDX = 2
  integer, parameter :: Nmax  = 4
    !! maximum #triangles associated with a node

contains

  module subroutine node_contain_pnt(this, pnt, res, with_edge)
    !! check if pnt lies within bnds of node
    class(node), intent(in)  :: this
    real,        intent(in)  :: pnt(2)
    logical,     intent(out) :: res
    logical,     intent(in), optional :: with_edge
      !! return res=true if pnt lies on a edge of node
      !! default: false

    integer :: i
    logical :: with_edge_, wprint_
    real, parameter :: eps = 1e-14

    with_edge_ = .false.
    if (present(with_edge)) with_edge_ = with_edge

    if (with_edge_) then
      if (pnt(X_IDX) >= this%bnds(X_IDX, LOWER) - eps .and. pnt(X_IDX) <= this%bnds(X_IDX, UPPER) + eps .and. &
        & pnt(Y_IDX) >= this%bnds(Y_IDX, LOWER) - eps .and. pnt(Y_IDX) <= this%bnds(Y_IDX, UPPER) + eps) then
        res =.true.
      else
        res = .false.
      end if
    else
      if (pnt(X_IDX) > this%bnds(X_IDX, LOWER) - eps .and. pnt(X_IDX) < this%bnds(X_IDX, UPPER) + eps .and. &
        & pnt(Y_IDX) > this%bnds(Y_IDX, LOWER) - eps .and. pnt(Y_IDX) < this%bnds(Y_IDX, UPPER) + eps) then
        res =.true.
      else
        res = .false.
      end if
    end if
  end subroutine

  module subroutine quadtree_init(this, g, Ntri_max, Nnodes)
    !! int quadtree
    class(quadtree),   intent(out) :: this
    type(triang_grid), intent(in)  :: g
    integer,           intent(in)  :: Ntri_max
      !! maximum #triangles associated with a node
    integer,           intent(in)  :: Nnodes
      !! maximum #nodes in quadtree

    integer    :: i, Ntri
    type(node) :: n1

    ! store grid data
    this%cell2vert = g%cell2vert
    this%vert      = g%vert
    Ntri           = size(g%cell2vert(1,:))
    call this%itr_vec%init(Ntri, c=5*Ntri, x=[(i, i=1, Ntri)])

    ! init quadtree
    this%Ntri_max = Ntri_max
    call this%nodes%init(0, c=Nnodes)

    ! init first node
    do i=1, 2
      n1%bnds(i,:) = [minval(g%vert(i,:)), maxval(g%vert(i,:))]
    end do
    n1%ichild = 0
    n1%ilower = 1
    n1%iupper = Ntri
    call this%nodes%push(n1)

    ! create quadtree
    do i=1, Nnodes
      if (i > this%nodes%n) exit
      call this%subdivide(this%nodes%d(i))
      if (this%nodes%n > Nnodes) then
        print *, "this%nodes%n", this%nodes%n
        call program_error('The given upper boundary for #nodes has been exceeded. Choose a larger Nnodes.')
      end if
    end do
  end subroutine

  module subroutine quadtree_subdivide(this, n)
    !! subdivide node into its 4 children if necessary

    class(quadtree), intent(inout) :: this
    type(node),      intent(inout) :: n
      !! node (parent)

    integer          :: i, ix, iy
    real             :: dr(2)
    type(node)       :: nchild(2,2)
    type(vector_int) :: itr_vec_

    if (n%iupper-n%ilower+1 <= this%Ntri_max) then
      !! leaf found
      n%ichild = 0
      return
    end if

    ! create children
    n%ichild = this%nodes%n + 1

    dr = 0.5*(n%bnds(:,2)-n%bnds(:,1))
    do ix=1, 2; do iy=1, 2
      nchild(ix,iy)%bnds(:,1) = n%bnds(:,1) + [real(ix)-1.0, real(iy)-1.0]*dr ! lower, left corner
      nchild(ix,iy)%bnds(:,2) = nchild(ix,iy)%bnds(:,1) + dr

      ! determine which triangles are within bnds of nchild(ix,iy)
      call itr_vec_%init(0, c=n%iupper-n%ilower+1)
      do i=n%ilower, n%iupper
        if (this%overlap(nchild(ix,iy), this%itr_vec%d(i))) call itr_vec_%push(this%itr_vec%d(i))
          !! push triangles overlapping with bnds of nchild
      end do

      ! store indices of overlapping triangles in this%itr_vec
      nchild(ix,iy)%ilower = this%itr_vec%n + 1
      nchild(ix,iy)%iupper = this%itr_vec%n + itr_vec_%n
      call this%itr_vec%push(itr_vec_%to_array())

      call this%nodes%push(nchild(ix,iy))
      call itr_vec_%destruct()
    end do; end do
  end subroutine

  module function quadtree_overlap(this, n, itriang) result(res)
    class(quadtree), intent(in) :: this
    type(node),      intent(in) :: n
      !! node
    integer,         intent(in) :: itriang
      !! idx of triangle
    logical                     :: res

    integer :: i, ivert
    real    :: vert(2,3), c(2)
    logical :: overlap

    ! check if any of the triangle's vertices are within the node
    res = .true.
    do i=1, 3
      ivert = this%cell2vert(i, itriang)
      call n%contain_pnt(this%vert(:,ivert), overlap)
      if (overlap) return
    end do

    ! get triangle's vertices
    do i=1, 3
      ivert     = this%cell2vert(i, itriang)
      vert(:,i) = this%vert(:,ivert)
    end do

    ! check if triangle's center is within node
    c = [sum(vert(X_IDX, :)), sum(vert(Y_IDX, :))]/3
    call n%contain_pnt(c, overlap)
    if (overlap) return

    ! check if there is an intersection between triangle and node. Check all possible line intersections.
    block
      integer :: dir, dir2, r, di(2), dj(2)
      real    :: dr(2), l1(2,2), l2(2,2)

      dr = n%bnds(:,2)-n%bnds(:,1)
      do i=1, 3
        ! get triangle line
        l1(:,1) = vert(:,i)
        if (i==3) then
          l1(:,2) = vert(:,1)
        else
          l1(:,2) = vert(:,i+1)
        end if

        ! get line through vertices of node
        do dir = 1, 2
          dir2 = 3 - dir

          di(dir ) = 1
          di(dir2) = 0
          dj(dir ) = 0
          dj(dir2) = 1

          do r = 0, 1
            l2(:,1) = n%bnds(:,1) + (r*dj)*dr
            l2(:,2) = l2(:,1)     + di*dr
            if (lines_intersect(l1,l2)) return
          end do
        end do
      end do

    end block

    ! check if any of the node's vertices are within the triangle
    block
      integer :: ix, iy
      real    :: dr(2), p(2)

      dr = n%bnds(:,2)-n%bnds(:,1)
      do ix=1, 2; do iy=1, 2
        p = n%bnds(:,1) + [ix-1, iy-1]*dr ! node's vertices
        if (tri_contain_pnt(vert, p)) return
      end do; end do
    end block

    res = .false.
  end function

  module function quadtree_lookup_pnt(this, pnt) result(icell)
    !! find idx of the triangle that contains the point pnt

    class(quadtree), intent(in) :: this
    real,            intent(in) :: pnt(2)
      !! point
    integer                     :: icell
      !! result: idx of triangle

    integer :: i, j, k, l, itr
    logical :: within_node
    real    :: vert(2,3)

    icell = 0 ! not found

    i = 1; k = 0
    do while(i <= this%nodes%n)
      call this%nodes%d(i)%contain_pnt(pnt, within_node, with_edge=.true.)
      if (within_node) then
        if (this%nodes%d(i)%ichild == 0) then
          !! leaf found

          ! iterate through triangles of node and find the correct one
          do j=this%nodes%d(i)%ilower, this%nodes%d(i)%iupper
            itr = this%itr_vec%d(j)
            do l=1, 3
              vert(:,l) = this%vert(:,this%cell2vert(l, itr))
            end do
            if (tri_contain_pnt(vert, pnt)) then
              !! triangle found
              icell = itr
              return
            end if
          end do
        else
          k = 0
          i = this%nodes%d(i)%ichild
        end if
      else
        i = i + 1
        k = k + 1
        if (k > 4) return ! a node cannot have more than 4 childern
      end if
    end do
  end function

  module subroutine quadtree_print(this)
    class(quadtree), intent(in) :: this

    integer :: i

    print '(A, I10, /)', 'total #nodes', this%nodes%n
    print '(A)', 'all leafs:'
    do i=1, this%nodes%n
      if (this%nodes%d(i)%ichild == 0) then
        print '(F6.2, x, F6.2, x, F6.2, x, F6.2, I3)', this%nodes%d(i)%bnds(1,1), this%nodes%d(i)%bnds(1,2), &
                                                      &this%nodes%d(i)%bnds(2,1), this%nodes%d(i)%bnds(2,2), &
                                                      &(this%nodes%d(i)%iupper - this%nodes%d(i)%ilower + 1)
      end if
    end do
  end subroutine

  module subroutine quadtree_output(this, fname)
    class(quadtree), intent(in) :: this
    character(*),    intent(in) :: fname
      !! file name

    integer :: iounit, ios, i

    open (newunit=iounit, file=fname, iostat=ios, action="WRITE")
    if (ios /= 0) call program_error("Error opening file")

    do i=1, this%nodes%n
      if (this%nodes%d(i)%ichild == 0) then
        write (iounit, '(ES25.16E3, x, ES25.16E3, x, ES25.16E3, x, ES25.16E3)'),      &
              & minval(this%nodes%d(i)%bnds(1,:)), minval(this%nodes%d(i)%bnds(2,:)), &
              & abs(this%nodes%d(i)%bnds(1,2)-this%nodes%d(i)%bnds(1,1)), abs(this%nodes%d(i)%bnds(2,2)-this%nodes%d(i)%bnds(2,1))
      end if
    end do

    close (iounit)
  end subroutine

  function tri_contain_pnt(vert, p) result(res)
    real, intent(in) :: vert(2,3)
      !! vertices of triangle
    real, intent(in) :: p(2)
      !! point
    logical          :: res

    real :: s, t, A

    res = .false.

    ! area of triangle
    A = 0.5 *(-vert(Y_IDX,2)*vert(X_IDX,3) + vert(Y_IDX,1)*(-vert(X_IDX,2)+vert(X_IDX,3)) + vert(X_IDX,1)*(vert(Y_IDX,2)-vert(Y_IDX,3)) + vert(X_IDX,2)*vert(Y_IDX,3))

    ! use barycentric coordinates to determine if p(:) is in triangle
    s = (vert(Y_IDX,1)*vert(X_IDX,3) - vert(X_IDX,1)*vert(Y_IDX,3) + (vert(Y_IDX,3)-vert(Y_IDX,1))*p(X_IDX) + (vert(X_IDX,1)-vert(X_IDX,3))*p(Y_IDX))/(2*A)
    t = (vert(X_IDX,1)*vert(Y_IDX,2) - vert(Y_IDX,1)*vert(X_IDX,2) + (vert(Y_IDX,1)-vert(Y_IDX,2))*p(X_IDX) + (vert(X_IDX,2)-vert(X_IDX,1))*p(Y_IDX))/(2*A)
    if (s > 0 .and. t > 0 .and. s+t < 1) res = .true.
  end function

  function lines_intersect(l1, l2) result(res)
    !! do the lines l1 and l2 intersect?

    real, intent(in) :: l1(2,2), l2(2,2)
      !! lines. first index: x, y idx. second index: 1st or 2nd pnt of line
    logical          :: res

    real :: s, t, a(2), b(2), c(2)

    a = l1(:,2)-l1(:,1)
    b = l2(:,2)-l2(:,1)
    c = l2(:,1)-l1(:,1)

    res = .false.
    if (cross_product_2d(a,b) == 0) return! lines parallel

    ! l1 = l1(:,1) + s*(l1(:,2)-l1(:,1)),  l2 = l2(:,1) + t*(l2(:,2)-l2(:,1))
    s = cross_product_2d(c, b)  * cross_product_2d(a, b) / abs(cross_product_2d(a, b))**2
    t = cross_product_2d(-c, a) * cross_product_2d(b, a) / abs(cross_product_2d(b, a))**2

    if (s>0 .and. s<1 .and. t>0 .and. t<1) res = .true.
  end function
end submodule
