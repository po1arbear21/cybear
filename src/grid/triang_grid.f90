#include "../util/macro.f90.inc"

module triang_grid_m

  use error_m
  use grid_m,   only: grid, IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL
  use math_m,   only: cross_product_2d
  use vector_m, only: vector_int

  implicit none

  private
  public triang_grid

  type, extends(grid) :: triang_grid
    !! unstructured 2d triangle grid
    real, allocatable :: vert(:,:)
      !! vertices:  [x_i, y_i] = vert(1:2,i).
      !! size: (2,nvert)
    integer, allocatable :: cell2vert(:,:)
      !! indices to cell's vertices: cell_i = [vert_1, vert_2, vert_3] = cell2vert(1:3,i)
      !! size(3,ncell)

    integer, allocatable :: edge2vert(:,:)
      !! size: (2,nedge)

    integer, allocatable :: vert2edge_n(:), vert2edge_i(:)
      !! input:  vertex index i
      !! output: edges' number for vertex i.
      !!         how many edges for vertex i? vert2edge_n(i+1)-vert2edge_n(i)
      !!         edges of vertex i? vert2edge_i(j), j=vert2edge_n(i),..,vert2edge_n(i+1)-1
      !! sizes: vert2edge_n(1:nvert+1), this%vert2edge_i(1:2*nedge))

    integer, allocatable :: vert2cell_n(:), vert2cell_i(:)
      !! input:  vertex index i
      !! output: cells' number for vertex i.
      !!         how many cells for vertex i? vert2cell_n(i+1)-vert2cell_n(i)
      !!         cells of vertex i? vert2cell_i(j), j=vert2cell_n(i),..,vert2cell_n(i+1)-1
      !! sizes: vert2cell_n(1:nvert+1), this%vert2cell_i(1:3*ncell))

    integer :: max_neighb(4,4)
      !! maximal neighbours wrt idx_types
      !! indices: (idx1_type, idx2_type), each idx_type = IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL
      !! max_neighb(idx1_type, idx2_type): idx1_type has how many neighbours of type idx2_type

  contains
    procedure :: init           => triang_grid_init
    procedure :: get_idx_bnd    => triang_grid_get_idx_bnd
    procedure :: get_vertex     => triang_grid_get_vertex
    procedure :: get_edge       => triang_grid_get_edge
    procedure :: get_face       => triang_grid_get_face
    procedure :: get_cell       => triang_grid_get_cell
    procedure :: get_surf       => triang_grid_get_surf
    procedure :: get_vol        => triang_grid_get_vol
    procedure :: get_max_neighb => triang_grid_get_max_neighb
    procedure :: get_neighb     => triang_grid_get_neighb
  end type

contains

  subroutine triang_grid_init(this, vert, icell)
    !! initialize 1D grid
    class(triang_grid), intent(out) :: this
    real,               intent(in)  :: vert(:,:)
      !! vertices:  [x_i, y_i] = vert(1:2,i).
      !! size: (dim=2,nvert)
    integer,            intent(in)  :: icell(:,:)
      !! indices to cell's vertices: cell_i = [vert_1, vert_2, vert_3] = icell(1:3,i)
      !! size(3,ncell)

    integer                       :: ncell, nvert, nedge, iv
    type(vector_int), allocatable :: edge_oneway(:)
      !! edge (iv<->iv2) is only saved at one vertex, either iv or iv2.
    type(vector_int), allocatable :: vert2edge_vec(:), vert2cell_vec(:)

    ASSERT(2 == size(vert,  dim=1))
    ASSERT(3 == size(icell, dim=1))

    ! init base
    call this%grid_init(2, 1, [2], 3)

    ncell = size(icell, dim=2)
    nvert = size(vert, dim=2)

    ! save grid data
    this%vert = vert
    this%cell2vert = icell

    ! create tmp edges: oneway only
    ! create tmp: vert2cell_vec
    allocate (edge_oneway(nvert), vert2cell_vec(nvert))
    do iv = 1, nvert
      call edge_oneway(iv)%init(0, c=8)
      call vert2cell_vec(iv)%init(0, c=8)
    end do

    block
      integer :: ic, iv_, iv2

      do ic = 1, ncell
        do iv_ = 1, 3
          iv  = icell(iv_,ic)
          iv2 = icell(merge(iv_+1, 1, (iv_<3)),ic)

          ! create edge iv <-> iv2
          if (.not. (edge_exists(edge_oneway(iv), iv2) .or. edge_exists(edge_oneway(iv2), iv))) then
            call edge_oneway(iv)%push(iv2)
          end if

          call vert2cell_vec(iv)%push(ic)
        end do
      end do
    end block

    ! count edges
    nedge = 0
    do iv = 1, nvert
      nedge = nedge + edge_oneway(iv)%n
    end do

    allocate (vert2edge_vec(nvert))
    do iv = 1, nvert
      call vert2edge_vec(iv)%init(0, c=8)
    end do

    ! create edge2vert
    ! create tmp vert2edge_vec
    block
      integer :: ie, iv2_, iv2

      allocate (this%edge2vert(2,nedge))
      ie=1
      do iv = 1, nvert
        do iv2_ = 1, edge_oneway(iv)%n
          iv2 = edge_oneway(iv)%d(iv2_)
          this%edge2vert(:,ie) = [iv, iv2]
          call vert2edge_vec(iv )%push(ie)
          call vert2edge_vec(iv2)%push(ie)
          ie = ie + 1
        end do
      end do
    end block

    ! create vert2edge
    allocate (this%vert2edge_n(nvert+1), this%vert2edge_i(2*nedge))
    this%vert2edge_n(1) = 1

    do iv = 1, nvert
      this%vert2edge_n(iv+1) = this%vert2edge_n(iv) + vert2edge_vec(iv)%n
      this%vert2edge_i(this%vert2edge_n(iv):this%vert2edge_n(iv+1)-1) = vert2edge_vec(iv)%to_array()
    end do

    ! create vert2cell
    allocate (this%vert2cell_n(nvert+1), this%vert2cell_i(3*ncell))
    this%vert2cell_n(1) = 1

    do iv = 1, nvert
      this%vert2cell_n(iv+1) = this%vert2cell_n(iv) + vert2cell_vec(iv)%n
      this%vert2cell_i(this%vert2cell_n(iv):this%vert2cell_n(iv+1)-1) = vert2cell_vec(iv)%to_array()
    end do

    ! compute maximal neighbours
    block
      integer :: n, ntmp, iv_, ie

      this%max_neighb(IDX_VERTEX,IDX_VERTEX) = maxval(this%vert2edge_n(2:nvert+1) - this%vert2edge_n(1:nvert)) +1
      this%max_neighb(IDX_VERTEX,IDX_EDGE  ) = this%max_neighb(IDX_VERTEX,IDX_VERTEX) -1
      this%max_neighb(IDX_VERTEX,IDX_CELL  ) = maxval(this%vert2cell_n(2:nvert+1) - this%vert2cell_n(1:nvert))

      this%max_neighb(IDX_EDGE,IDX_VERTEX) = 2 ! = sum(this%face_dim)
      this%max_neighb(IDX_EDGE,IDX_CELL  ) = 2

      this%max_neighb(IDX_CELL,IDX_VERTEX) = 3 ! = this%cell_dim
      this%max_neighb(IDX_CELL,IDX_EDGE  ) = 3
      this%max_neighb(IDX_CELL,IDX_CELL  ) = 3

      ! find this%max_neighb(IDX_EDGE,IDX_EDGE)
      n = 0
      do ie = 1, nedge
        ntmp = 0
        do iv_ = 1, 2
          iv   = this%edge2vert(iv_,ie)
          ntmp = ntmp + this%vert2edge_n(iv+1)-this%vert2edge_n(iv)
        end do
        n = max(n, ntmp-1)
      end do
      this%max_neighb(IDX_EDGE,IDX_EDGE) = n

      ! 2D: face=edge
      this%max_neighb(IDX_FACE,IDX_FACE) = this%max_neighb(IDX_EDGE,IDX_EDGE)
      this%max_neighb(:,IDX_FACE) = this%max_neighb(:,IDX_EDGE)
      this%max_neighb(IDX_FACE,:) = this%max_neighb(IDX_EDGE,:)
    end block

  contains
    logical function edge_exists(edge, iv0)
      type(vector_int), intent(in) :: edge
      integer,          intent(in) :: iv0

      integer :: j

      edge_exists = (0 < count([(iv0 == edge%d(j), j=1,edge%n)]))
    end function
  end subroutine

  subroutine triang_grid_get_idx_bnd(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in)  :: idx_dir
      !! index direction for edges and faces. (must be 0 for IDX_VERTEX and IDX_CELL as always)
    integer,            intent(out) :: idx_bnd(:)
      !! output: upper bound for each index. size: (idx_dim=1)

    ASSERT(this%idx_allowed(idx_type, idx_dir))
    ASSERT(size(idx_bnd) == this%idx_dim)

    IGNORE(idx_dir)

    select case (idx_type)
      case (IDX_VERTEX)
        idx_bnd = size(this%vert,      dim=2)
      case (IDX_EDGE)
        idx_bnd = size(this%edge2vert, dim=2)
      case (IDX_FACE)
        idx_bnd = size(this%edge2vert, dim=2)
      case (IDX_CELL)
        idx_bnd = size(this%cell2vert, dim=2)
    end select
  end subroutine

  subroutine triang_grid_get_vertex(this, idx, p)
    !! get single vertex: from grid indices to coordinates
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! vertex' grid indices. size: (idx_dim=1)
    real,               intent(out) :: p(:)
      !! output: vertex' coordinates. size: (dim=2)

    ASSERT(this%idx_allowed(IDX_VERTEX, 0, idx=idx))
    ASSERT(size(p) == this%dim)

    p = this%vert(:,idx(1))
  end subroutine

  subroutine triang_grid_get_edge(this, idx, idx_dir, p)
    !! get single edge: from grid indices to coordinates
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! edge's indices. size: (1)
    integer,            intent(in)  :: idx_dir
      !! edge's index direction
    real,               intent(out) :: p(:,:)
      !! output: edge's coordinates. size: (dim=2, 2)

    integer :: iv

    ASSERT(this%idx_allowed(IDX_EDGE, idx_dir, idx=idx))
    ASSERT(all(shape(p) == [this%dim, 2]))

    IGNORE(idx_dir)

    do iv = 1, 2
      p(:,iv) = this%vert(:,this%edge2vert(iv,idx(1)))
    end do
  end subroutine

  subroutine triang_grid_get_face(this, idx, idx_dir, p)
    !! get single face: from grid indices to coordinates
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! face's indices. size: (idx_dim=1)
    integer,            intent(in)  :: idx_dir
      !! face's index direction
    real,               intent(out) :: p(:,:)
      !! output: face's coordinates. size: (dim=2, 2)

    ASSERT(this%idx_allowed(IDX_FACE, idx_dir, idx=idx))
    ASSERT(all(shape(p) == [this%dim, this%face_dim]))

    call this%get_edge(idx, idx_dir, p)
  end subroutine

  subroutine triang_grid_get_cell(this, idx, p)
    !! get single cell: from grid indices to coordinates
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! cell's indices. size: (idx_dim=1)
    real,               intent(out) :: p(:,:)
      !! output: cell's coordinates. size: (dim=2, cell_dim=3)

    integer :: iv

    ASSERT(this%idx_allowed(IDX_CELL, 0, idx=idx))
    ASSERT(all(shape(p) == [this%dim, this%cell_dim]))

    do iv = 1, 3
      p(:,iv) = this%vert(:,this%cell2vert(iv,idx(1)))
    end do
  end subroutine

  function triang_grid_get_surf(this, idx, idx_dir) result(surf)
    !! get single face's surface
    class(triang_grid), intent(in) :: this
    integer,            intent(in) :: idx(:)
      !! face's indices. size: (idx_dim=1)
    integer,            intent(in) :: idx_dir
      !! face's index direction
    real                           :: surf
      !! output: size of face

    real :: p(2,2)

    ASSERT(this%idx_allowed(IDX_FACE, idx_dir, idx=idx))

    call this%get_edge(idx, idx_dir, p)
    surf = norm2(p(:,1) - p(:,2))
  end function

  function triang_grid_get_vol(this, idx) result(vol)
    !! get single cell's volume
    class(triang_grid), intent(in) :: this
    integer,            intent(in) :: idx(:)
      !! cell's indices. size: (idx_dim=1)
    real                           :: vol
      !! return cell volume

    integer :: iv
    real    :: v(2,3)

    ASSERT(this%idx_allowed(IDX_CELL, 0, idx=idx))

    ! vertices
    do iv = 1, 3
      v(:,iv) = this%vert(:,this%cell2vert(iv,idx(1)))
    end do

    vol = 0.5 * abs(cross_product_2d(v(:,3)-v(:,1), v(:,2)-v(:,1)))
  end function

  function triang_grid_get_max_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir) result(max_neighb)
    !! get maximal number of nearest neighbours
    class(triang_grid), intent(in) :: this
    integer,            intent(in) :: idx1_type
      !! first index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in) :: idx1_dir
      !! first index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,            intent(in) :: idx2_type
      !! neighbour index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in) :: idx2_dir
      !! neighbour index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer                        :: max_neighb
      !! return maximal number of nearest neighbours

    ASSERT(this%idx_allowed(idx1_type, idx1_dir))
    ASSERT(this%idx_allowed(idx2_type, idx2_dir))

    IGNORE(this)
    IGNORE(idx1_dir)
    IGNORE(idx2_dir)

    max_neighb = this%max_neighb(idx1_type,idx2_type)
  end function

  subroutine triang_grid_get_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir, idx1, j, idx2, status)
      !! get j-th neighbour.
      !!
      !! j: we count neighbors from 1,2,...,N.
      !! N: depends on idx1 (e.g. boundary nodes might have fewer neighbors).
      !! status: indicates if j-th neighbor exists
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx1_type
      !! first index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in)  :: idx1_dir
      !! first index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,            intent(in)  :: idx2_type
      !! neighbour index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in)  :: idx2_dir
      !! neighbour index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,            intent(in)  :: idx1(:)
      !! first indices. size: (idx_dim)
    integer,            intent(in)  :: j
      !! j-th neighbor
    integer,            intent(out) :: idx2(:)
      !! output neighbour indices. size: (idx_dim)
    logical,            intent(out) :: status
      !! does j-th neighbor exist?

    integer :: max_neighb

    ASSERT(this%idx_allowed(idx1_type, idx1_dir, idx=idx1))
    ASSERT(this%idx_allowed(idx2_type, idx2_dir))
    ASSERT(size(idx2) == size(idx1))

    max_neighb = this%get_max_neighb(idx1_type, idx1_dir, idx2_type, idx2_dir)

    status = .false.
    if ((j < 1) .or. (j > max_neighb)) return

    if (idx1_type == IDX_VERTEX) then
      block
        integer :: iv1

        iv1 = idx1(1)

        if (idx2_type == IDX_VERTEX) then
          block
            integer :: ie, ie_, iv2, iv2_

            if (j <= this%vert2edge_n(iv1+1) - this%vert2edge_n(iv1)) then  ! j must be smaller/equal than number of edges
              ie_ = this%vert2edge_n(iv1) + j-1
              ie  = this%vert2edge_i(ie_)
              do iv2_ = 1, 2
                iv2 = this%edge2vert(iv2_,ie)
                if (iv2 == iv1) cycle
                idx2 = iv2
              end do
              status = .true.
            end if
          end block

        else if ((idx2_type == IDX_EDGE) .or. (idx2_type == IDX_FACE)) then
          block
            integer :: ie_

            if (j <= this%vert2edge_n(iv1+1) - this%vert2edge_n(iv1)) then  ! j must be smaller/equal than number of edges
              ie_  = this%vert2edge_n(iv1) + j-1
              idx2 = this%vert2edge_i(ie_)
            end if
          end block

        else if (idx2_type == IDX_CELL) then
          block
            integer :: ic_

            if (j <= this%vert2edge_n(iv1+1) - this%vert2edge_n(iv1)) then  ! j must be smaller/equal than number of cells
              ic_  = this%vert2cell_n(iv1) + j-1
              idx2 = this%vert2cell_i(ic_)
            end if
          end block
        end if            ! idx2_type
      end block           ! idx1_type==IDX_VERTEX

    else if ((idx1_type == IDX_EDGE) .or. (idx1_type == IDX_FACE)) then
      block
        integer :: ie1

        ie1 = idx1(1)

        if (idx2_type == IDX_VERTEX) then
          if (j < 3) idx2 = this%edge2vert(j,ie1)

        else if ((idx2_type == IDX_EDGE) .or. (idx2_type == IDX_FACE)) then
          block
            integer :: iv1, iv1_, ie2, ie2_, ne2

            ! how many edges we encountered this far
            ne2 = 0

            OUTER: do iv1_ = 1, 2
              iv1 = this%edge2vert(iv1_,ie1)
              ! j must be smaller/equal than number of edges for that specific vertex
              ! (minus 1, b.c. of self counting edge_1)
              if (j-ne2 <= this%vert2edge_n(iv1+1) - this%vert2edge_n(iv1)-1) then
                do ie2_ = this%vert2edge_n(iv1), this%vert2edge_n(iv1+1)-1
                  ie2 = this%vert2edge_i(ie2_)
                  if (ie1 == ie2) cycle
                  ne2 = ne2 + 1
                  if (ne2 == j) then
                    idx2   = ie2
                    status = .true.
                    exit OUTER
                  end if
                end do
              else
                ne2 = ne2 + this%vert2edge_n(iv1+1) - this%vert2edge_n(iv1)-1
              end if
            end do OUTER
          end block

        else if (idx2_type == IDX_CELL) then
          block
            integer :: iv1(2), ic2, ic2_, iv2_, nc2

            ! make sure iv1(1) has fewer neighb cells than iv1(2) as we only loop over one vertex' cells
            iv1 = this%edge2vert(:,ie1)
            if (this%vert2cell_n(iv1(2)+1)-this%vert2cell_n(iv1(2)) < this%vert2cell_n(iv1(1)+1) - this%vert2cell_n(iv1(1))) then
              iv1 = [iv1(2), iv1(1)]
            end if

            ! how many cells we encountered this far
            nc2 = 0

            do ic2_ = this%vert2cell_n(iv1(1)), this%vert2cell_n(iv1(1)+1)-1
              ic2 = this%vert2cell_i(ic2_)

              ! does iv1(1)'s neighb cell ic2 contain iv1(2) aka is edge ie1 part of ic2?
              if (count([(this%cell2vert(iv2_,ic2) == iv1(2), iv2_ = 1, 3)]) > 0) then
                nc2 = nc2 + 1
                if (nc2 == j) then
                  idx2   = ic2
                  status = .true.
                  exit
                end if
              end if
            end do
          end block
        end if
      end block

    else if (idx1_type == IDX_CELL) then
      block
        integer :: ic1

        ic1 = idx1(1)

        if (idx2_type == IDX_VERTEX) then
          idx2   = this%cell2vert(j,ic1)
          status = .true.

        else if ((idx2_type == IDX_EDGE) .or. (idx2_type == IDX_FACE)) then
          block
            integer :: iv1, iv1_, iv11, iv11_, ie2, ie2_

            iv1_  = j
            iv11_ = merge(iv1_+1, 1, iv1_<3)
            iv1   = this%cell2vert(iv1_, ic1)
            iv11  = this%cell2vert(iv11_,ic1)
            do ie2_ = this%vert2edge_n(iv1), this%vert2edge_n(iv1+1)-1
              ie2 = this%vert2edge_i(ie2_)
              if (any(iv11 == this%edge2vert(:,ie2))) then
                idx2   = ie2
                status = .true.
                exit
              end if
            end do
          end block

        else if (idx2_type == IDX_CELL) then
          block
            integer :: iv1, iv1_, iv11, iv11_, ic2, ic2_, nc2

            nc2 = 0
            do iv1_ = 1, 3
              iv11_ = merge(iv1_+1, 1, iv1_<3)
              iv1  = this%cell2vert(iv1_, ic1)
              iv11 = this%cell2vert(iv11_,ic1)
              do ic2_ = this%vert2cell_n(iv1), this%vert2cell_n(iv1+1)-1
                ic2 = this%vert2cell_i(ic2_)
                if (ic2 == ic1) cycle
                if (any(iv11 == this%cell2vert(:, ic2))) then
                  nc2 = nc2 + 1
                  if (nc2 == j) then
                    idx2   = ic2
                    status = .true.
                    exit
                  end if
                end if
              end do
            end do
          end block
        end if            ! idx2_type
      end block
    end if                ! idx1_type
  end subroutine

end module
