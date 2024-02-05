m4_include(../util/macro.f90.inc)

module triang_grid_m

  use error_m,         only: assert_failed, program_error
  use grid_m,          only: grid, IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL
  use json_m,          only: json_object
  use math_m,          only: cross_product_2d
  use normalization_m, only: denorm, norm
  use output_file_m,   only: output_file
  use plotmtv_m,       only: plotmtv, plotset_options
  use qsort_m,         only: qsort
  use vector_m,        only: vector_int, vector_real

  implicit none

  private
  public triang_grid

  type, extends(grid) :: triang_grid
    !! unstructured 2D triangle grid

    integer :: nvert
      !! number of vertices
    integer :: ncell
      !! number of cells
    integer :: nedge
      !! number of edges

    real, allocatable :: vert(:,:)
      !! vertex coordinates: [x_i, y_i] = vert(1:2,i); size = (2, nvert)

    integer, allocatable :: vert2edge_ia(:), vert2edge_ja(:)
      !! get edge neighbours of vertex (number varies, use CSR format)
    integer, allocatable :: vert2cell_ia(:), vert2cell_ja(:)
      !! get cell neighbours of vertex (number varies, use CSR format)
    integer, allocatable :: edge2vert(:,:)
      !! get vertex neighbours of edge (endpoints); size = (2, nedge)
    integer, allocatable :: edge2edge_ia(:), edge2edge_ja(:)
      !! get edge neighbours of edge (number varies, use CSR format)
    integer, allocatable :: edge2cell(:,:)
      !! get cell neighbours of edge (2 triangles share 1 edge); size = (2, nedge)
    integer, allocatable :: cell2vert(:,:)
      !! get vertex neighbours of cell (triangle corners); size = (3, ncell)
    integer, allocatable :: cell2edge(:,:)
      !! get edge neighbours of cell (triangle sides); size = (3, ncell)
    integer, allocatable :: cell2cell(:,:)
      !! get cell neighbours of cell; size = (3, ncell)

    integer :: max_neighb(4,4)
      !! maximal number of neighbours between two index types (IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL)

    type(quadtree), allocatable :: qtree
      !! quadtree to efficiently search cells based on point coordinates
  contains
    procedure :: init           => triang_grid_init
    procedure :: get_icell      => triang_grid_get_icell
    procedure :: get_idx_bnd_n  => triang_grid_get_idx_bnd_n
    procedure ::                   triang_grid_get_idx_bnd_1
    generic   :: get_idx_bnd    => triang_grid_get_idx_bnd_1
    procedure :: get_vertex     => triang_grid_get_vertex
    procedure :: get_edge       => triang_grid_get_edge
    procedure :: get_face       => triang_grid_get_face
    procedure :: get_cell       => triang_grid_get_cell
    procedure :: get_len        => triang_grid_get_len
    procedure :: get_surf       => triang_grid_get_surf
    procedure :: get_vol        => triang_grid_get_vol
    procedure :: get_max_neighb => triang_grid_get_max_neighb
    procedure :: get_neighb     => triang_grid_get_neighb
    procedure :: get_adjoint    => triang_grid_get_adjoint
    procedure :: output         => triang_grid_output
    procedure :: output_plotmtv => triang_grid_output_plotmtv
    procedure :: output_matlab  => triang_grid_output_matlab
    procedure :: read_matlab    => triang_grid_read_matlab
  end type

  type node
    integer :: ilower, iupper
      !! vars give indices of all triangles that are within the node,
      !!  i.e. quadtree%itr_vec%d(ilower:iupper) == all triangle indices
    integer :: ichild
      !! index of childrens in quadtree%n(:)
      !! ichild = 0 (no children)
      !!        = i (nodes i, i+1, i+2, i+3 are children)
    real    :: bnds(2,2)
      !! axis-aligned boundaries
      !! 1st index: x-  (=: 1) or y-range   (=: 2)
      !! 2nd index: min (=: 1) or max value (=: 2)
  contains
    procedure :: contain_pnt => node_contain_pnt
  end type

  m4_define({T},{node})
  m4_include(../util/vector_def.f90.inc)

  type quadtree
    integer, allocatable :: itr(:)
      !! indices of all triangles
    integer, allocatable :: cell2vert(:,:)
      !! indices to triangle's vertices: cell_i = [vert_1, vert_2, vert_3] = cell2vert(1:3,i)
      !! size(3,ncell)
    real,    allocatable :: vert(:,:)
      !! vertices:  [x_i, y_i] = vert(1:2,i).
    integer              :: Ntri_max
      !! maximum #triangles associated with a node. default: 4
    type(vector_node)    :: nodes
      !! nodes of quadtree
    type(vector_int)     :: itr_vec
      !! vector containing the triangle indices of all nodes
  contains
    procedure :: init       => quadtree_init
    procedure :: subdivide  => quadtree_subdivide
    procedure :: overlap    => quadtree_overlap
    procedure :: lookup_pnt => quadtree_lookup_pnt
    procedure :: print      => quadtree_print
  end type

  ! interface to submodule routines in grid/quadtree.f90
  interface
    module subroutine node_contain_pnt(this, pnt, res, with_edge)
      !! check if pnt lies within bnds of node
      class(node), intent(in)  :: this
      real,        intent(in)  :: pnt(2)
      logical,     intent(out) :: res
      logical,     intent(in), optional :: with_edge
    end subroutine

    module subroutine quadtree_init(this, g, Ntri_max, Nnodes)
      !! int quadtree
      class(quadtree),   intent(out) :: this
      type(triang_grid), intent(in)  :: g
      integer,           intent(in)  :: Ntri_max
        !! maximum #triangles associated with a node
      integer,           intent(in)  :: Nnodes
        !! maximum #nodes in quadtree
    end subroutine

    module subroutine quadtree_subdivide(this, n)
      !! subdivide node into its 4 children

      class(quadtree), intent(inout) :: this
      type(node),      intent(inout) :: n
    end subroutine

    module function quadtree_overlap(this, n, itriang) result(res)
      class(quadtree), intent(in) :: this
      type(node),      intent(in) :: n
        !! node
      integer,         intent(in) :: itriang
        !! idx of triangle
      logical                     :: res
    end function

    module function quadtree_lookup_pnt(this, pnt) result(icell)
      !! find idx of the triangle that contains the point pnt

      class(quadtree), intent(in) :: this
      real,            intent(in) :: pnt(2)
        !! point
      integer                     :: icell
        !! result: idx of triangle
    end function

    module subroutine quadtree_print(this)
      class(quadtree), intent(in) :: this
    end subroutine
  end interface

  integer, parameter :: INVALID_NEIGHBOUR = huge(1)
    !! mark neighbour as invalid (huge positive number so invalid indices get sorted to the back)

contains

  m4_define({T},{node})
  m4_include(../util/vector_imp.f90.inc)

  subroutine triang_grid_init(this, name, vert, cells)
    !! initialize triangle grid
    class(triang_grid), intent(out) :: this
    character(*),       intent(in)  :: name
      !! grid name
    real,               intent(in)  :: vert(:,:)
      !! x,y vertex coordinates; size = (2, nvert)
    integer,            intent(in)  :: cells(:,:)
      !! 3 vertex indices per cell; size = (3, ncell)

    integer                       :: j, k, icell, iedge1, iedge2, ivert1, ivert2, nmax
    integer,          allocatable :: edges(:,:)
    type(vector_int)              :: edge2vert, edge2cell, edge2edge
    type(vector_int), allocatable :: vert2edge(:), vert2cell(:)

    ! init base
    call this%grid_init(name, 2, 1, [2], 3, [3])

    ! save vertices and cell-vertex table
    m4_assert(size(vert,  dim = 1) == 2)
    m4_assert(size(cells, dim = 1) == 3)
    this%nvert     = size(vert,  dim = 2)
    this%ncell     = size(cells, dim = 2)
    this%vert      = vert
    this%cell2vert = cells

    ! allocate memory
    allocate (this%cell2edge(3,this%ncell), this%cell2cell(3,this%ncell), edges(this%nvert,this%nvert), source = INVALID_NEIGHBOUR)
    allocate (vert2edge(this%nvert), vert2cell(this%nvert))
    call edge2vert%init(0, c = 4 * this%nvert)
    call edge2cell%init(0, c = 2 * this%nvert)
    do ivert1 = 1, this%nvert
      call vert2edge(ivert1)%init(0, c = 6)
      call vert2cell(ivert1)%init(0, c = 6)
    end do

    ! create neighbour tables
    this%nedge = 0
    do icell = 1, this%ncell
      do j = 1, 3
        ! update vertex-cell table, use j as vertex index
        call vert2cell(cells(j,icell))%push(icell)

        ! get endpoints of edge, use j as edge index
        k = mod(j, 3) + 1
        ivert1 = min(cells(j,icell), cells(k,icell))
        ivert2 = max(cells(j,icell), cells(k,icell))

        ! check whether edge exists already
        if (edges(ivert1,ivert2) /= INVALID_NEIGHBOUR) then
          edge2cell%d(2 * edges(ivert1,ivert2)) = icell ! fill in placeholder
          this%cell2edge(j,icell) = edges(ivert1,ivert2)
        else
          ! add new edge
          this%nedge = this%nedge + 1
          edges(ivert1,ivert2) = this%nedge

          ! update tables
          this%cell2edge(j,icell) = this%nedge
          call edge2vert%push(ivert1)
          call edge2vert%push(ivert2)
          call edge2cell%push(icell)
          call edge2cell%push(INVALID_NEIGHBOUR) ! placeholder for possible second cell
          call vert2edge(ivert1)%push(this%nedge)
          call vert2edge(ivert2)%push(this%nedge)
        end if
      end do
    end do
    this%edge2vert = reshape(edge2vert%d(1:edge2vert%n), [2, this%nedge])
    this%edge2cell = reshape(edge2cell%d(1:edge2cell%n), [2, this%nedge])

    ! compress vertex-cell, vertex-edge neighbour tables (CSR format)
    allocate (this%vert2edge_ia(this%nvert+1), this%vert2edge_ja(2*this%nedge))
    allocate (this%vert2cell_ia(this%nvert+1), this%vert2cell_ja(3*this%ncell))
    associate (eia => this%vert2edge_ia, eja => this%vert2edge_ja, cia => this%vert2cell_ia, cja => this%vert2cell_ja)
      eia(1) = 1
      cia(1) = 1
      do ivert1 = 1, this%nvert
        eia(ivert1+1) = eia(ivert1) + vert2edge(ivert1)%n
        cia(ivert1+1) = cia(ivert1) + vert2cell(ivert1)%n
        eja(eia(ivert1):eia(ivert1+1)-1) = vert2edge(ivert1)%d(1:vert2edge(ivert1)%n)
        cja(cia(ivert1):cia(ivert1+1)-1) = vert2cell(ivert1)%d(1:vert2cell(ivert1)%n)
      end do
    end associate

    ! edge-edge neighbour table (CSR format)
    allocate (this%edge2edge_ia(this%nedge+1))
    call edge2edge%init(0, c = 16 * this%nedge)
    this%edge2edge_ia(1) = 1
    nmax = 0
    do iedge1 = 1, this%nedge
      do j = 1, 2
        ivert1 = this%edge2vert(j,iedge1)
        do k = this%vert2edge_ia(ivert1), this%vert2edge_ia(ivert1 + 1) - 1
          iedge2 = this%vert2edge_ja(k)
          if (iedge2 == iedge1) cycle ! don't include edge itself in list of neighbours
          call edge2edge%push(iedge2)
        end do
      end do
      this%edge2edge_ia(iedge1+1) = edge2edge%n + 1
      nmax = max(nmax, edge2edge%n - this%edge2edge_ia(iedge1) + 1)
    end do
    this%edge2edge_ja = edge2edge%to_array()

    ! cell-cell neighbour table
    do icell = 1, this%ncell
      do j = 1, 3
        iedge1 = this%cell2edge(j,icell)
        this%cell2cell(j,icell) = this%edge2cell(1,iedge1)
        if (this%cell2cell(j,icell) == icell) this%cell2cell(j,icell) = this%edge2cell(2,iedge1)
      end do
    end do

    ! sort neighbour tables (except cell2vert)
    do ivert1 = 1, this%nvert
      call qsort(this%vert2edge_ja(this%vert2edge_ia(ivert1):this%vert2edge_ia(ivert1+1)-1))
      call qsort(this%vert2cell_ja(this%vert2cell_ia(ivert1):this%vert2cell_ia(ivert1+1)-1))
    end do
    do iedge1 = 1, this%nedge
      call qsort(this%edge2edge_ja(this%edge2edge_ia(iedge1):this%edge2edge_ia(iedge1+1)-1))
      call qsort(this%edge2cell(:,iedge1))
    end do
    do icell = 1, this%ncell
      call qsort(this%cell2edge(:,icell))
      call qsort(this%cell2cell(:,icell))
    end do

    ! get maximum number of neighbours
    this%max_neighb(IDX_VERTEX,IDX_VERTEX) = maxval(this%vert2edge_ia(2:this%nvert+1) - this%vert2edge_ia(1:this%nvert))
    this%max_neighb(IDX_VERTEX,IDX_EDGE  ) = this%max_neighb(IDX_VERTEX,IDX_VERTEX)
    this%max_neighb(IDX_VERTEX,IDX_CELL  ) = maxval(this%vert2cell_ia(2:this%nvert+1) - this%vert2cell_ia(1:this%nvert))
    this%max_neighb(IDX_EDGE,  IDX_VERTEX) = 2
    this%max_neighb(IDX_EDGE,  IDX_EDGE  ) = 0
    this%max_neighb(IDX_EDGE,  IDX_EDGE  ) = nmax
    this%max_neighb(IDX_EDGE,  IDX_CELL  ) = 2
    this%max_neighb(IDX_CELL,           :) = 3
    this%max_neighb(       :,  IDX_FACE  ) = this%max_neighb(:,IDX_EDGE)
    this%max_neighb(IDX_FACE,           :) = this%max_neighb(IDX_EDGE,:)
  end subroutine

  subroutine triang_grid_get_icell(this, pnt, icell, Nnodes, Ntri_max)
    !! get cell idx of triangle that contains point pnt
    class(triang_grid), intent(inout) :: this
    real,               intent(in)    :: pnt(2)
      !! pnt = [x, y]
    integer,            intent(out)   :: icell
      !! idx of triangle
    integer, optional,  intent(in)    :: Ntri_max
      !! maximum #triangles associated with a node. default: 4
    integer, optional,  intent(in)    :: Nnodes
      !! maximum #nodes in quadtree. default: 256

    integer :: Nnodes_, Ntri_max_

    if (.not. allocated(this%qtree)) then
      Ntri_max_ = 4
      Nnodes_   = 256
      if (present(Ntri_max)) Ntri_max_ = Ntri_max
      if (present(Nnodes  )) Nnodes_   = Nnodes

      allocate (this%qtree)
      call this%qtree%init(this, Nnodes_, Ntri_max_)
    end if

    icell = this%qtree%lookup_pnt(pnt)
    if (icell == 0) call program_error("Unable to find cell")
  end subroutine

  subroutine triang_grid_get_idx_bnd_n(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in)  :: idx_dir
      !! index direction for edges and faces. (must be 0 for IDX_VERTEX and IDX_CELL as always)
    integer,            intent(out) :: idx_bnd(:,:)
      !! output: lower/upper bound for each index. size: (2, idx_dim = 1)

    m4_assert(size(idx_bnd,1) == 2)
    m4_assert(size(idx_bnd,2) == this%idx_dim)

    call this%get_idx_bnd(idx_type, idx_dir, idx_bnd(:,1))
  end subroutine

  subroutine triang_grid_get_idx_bnd_1(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in)  :: idx_dir
      !! index direction for edges and faces. (must be 0 for IDX_VERTEX and IDX_CELL as always)
    integer,            intent(out) :: idx_bnd(:)
      !! output: lower/upper bound for each index. (2)

    m4_assert(size(idx_bnd) == 2)
    m4_assert(this%idx_allowed(idx_type, idx_dir))

    idx_bnd(1) = 1
    select case (idx_type)
      case (IDX_VERTEX)
        idx_bnd(2) = this%nvert
      case (IDX_EDGE)
        idx_bnd(2) = this%nedge
      case (IDX_FACE)
        idx_bnd(2) = this%nedge
      case (IDX_CELL)
        idx_bnd(2) = this%ncell
    end select
  end subroutine

  subroutine triang_grid_get_vertex(this, idx, p)
    !! get single vertex: from grid indices to coordinates
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! vertex grid indices; size = (idx_dim=1)
    real,               intent(out) :: p(:)
      !! output: vertex coordinates; size = (dim=2)

    m4_assert(this%idx_allowed(IDX_VERTEX, 0, idx=idx))
    m4_assert(size(p) == this%dim)

    p = this%vert(:,idx(1))
  end subroutine

  subroutine triang_grid_get_edge(this, idx, idx_dir, p)
    !! get single edge: from grid indices to coordinates
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! edge indices; size = (idx_dim=1)
    integer,            intent(in)  :: idx_dir
      !! edge index direction = 1
    real,               intent(out) :: p(:,:)
      !! output edge coordinates; size = (dim=2, 2)

    m4_assert(this%idx_allowed(IDX_EDGE, idx_dir, idx=idx))
    m4_assert(all(shape(p) == [this%dim, 2]))

    p(:,1) = this%vert(:,this%edge2vert(1,idx(1)))
    p(:,2) = this%vert(:,this%edge2vert(2,idx(1)))
  end subroutine

  subroutine triang_grid_get_face(this, idx, idx_dir, p)
    !! get single face: from grid indices to coordinates
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! face indices; size = (idx_dim=1)
    integer,            intent(in)  :: idx_dir
      !! face index direction = 1
    real,               intent(out) :: p(:,:)
      !! output face coordinates; size = (dim=2, 2)

    m4_assert(this%idx_allowed(IDX_FACE, idx_dir, idx=idx))
    m4_assert(all(shape(p) == [this%dim, this%face_nvert(1)]))

    call this%get_edge(idx, idx_dir, p)
  end subroutine

  subroutine triang_grid_get_cell(this, idx, p)
    !! get single cell: from grid indices to coordinates
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! cell indices; size = (idx_dim=1)
    real,               intent(out) :: p(:,:)
      !! output: cell coordinates; size = (dim=2, cell_nvert=3)

    m4_assert(this%idx_allowed(IDX_CELL, 0, idx=idx))
    m4_assert(all(shape(p) == [this%dim, this%cell_nvert]))

    p(:,1) = this%vert(:,this%cell2vert(1,idx(1)))
    p(:,2) = this%vert(:,this%cell2vert(2,idx(1)))
    p(:,3) = this%vert(:,this%cell2vert(3,idx(1)))
  end subroutine

  function triang_grid_get_len(this, idx, idx_dir) result(len)
    !! get edge length
    class(triang_grid), intent(in) :: this
    integer,            intent(in) :: idx(:)
      !! edge indices (idx_dim=1)
    integer,            intent(in) :: idx_dir
      !! edge direction = 1
    real                           :: len
      !! return edge length

    real :: p(2,2)

    m4_assert(this%idx_allowed(IDX_EDGE, idx_dir, idx=idx))

    call this%get_edge(idx, idx_dir, p)
    len = norm2(p(:,1) - p(:,2))
  end function

  function triang_grid_get_surf(this, idx, idx_dir) result(surf)
    !! get face area
    class(triang_grid), intent(in) :: this
    integer,            intent(in) :: idx(:)
      !! face indices; size = (idx_dim=1)
    integer,            intent(in) :: idx_dir
      !! face index direction = 1
    real                           :: surf
      !! return face area

    surf = this%get_len(idx, idx_dir)
  end function

  function triang_grid_get_vol(this, idx) result(vol)
    !! get cell volume
    class(triang_grid), intent(in) :: this
    integer,            intent(in) :: idx(:)
      !! cell indices; size = (idx_dim=1)
    real                           :: vol
      !! return cell volume

    real :: v(2,3)

    m4_assert(this%idx_allowed(IDX_CELL, 0, idx=idx))

    call this%get_cell(idx, v)
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

    m4_assert(this%idx_allowed(idx1_type, idx1_dir))
    m4_assert(this%idx_allowed(idx2_type, idx2_dir))

    m4_ignore(this)
    m4_ignore(idx1_dir)
    m4_ignore(idx2_dir)

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

    integer :: i1, i2, k

    m4_assert(this%idx_allowed(idx1_type, idx1_dir, idx=idx1))
    m4_assert(this%idx_allowed(idx2_type, idx2_dir))
    m4_assert(size(idx2) == this%idx_dim)

    status = .false.
    idx2   = 0
    if ((j < 1) .or. (j > this%max_neighb(idx1_type,idx2_type))) return

    i1 = idx1(1)

    select case (idx1_type)
    case (IDX_VERTEX)
      select case (idx2_type)
      case (IDX_VERTEX)
        if (j > this%vert2edge_ia(i1 + 1) - this%vert2edge_ia(i1)) return
        k  = this%vert2edge_ja(this%vert2edge_ia(i1) + j - 1)
        i2 = this%edge2vert(1,k)
        if (i2 == i1) i2 = this%edge2vert(2,k)
      case (IDX_EDGE, IDX_FACE)
        if (j > this%vert2edge_ia(i1 + 1) - this%vert2edge_ia(i1)) return
        i2 = this%vert2edge_ja(this%vert2edge_ia(i1) + j - 1)
      case (IDX_CELL)
        if (j > this%vert2cell_ia(i1 + 1) - this%vert2cell_ia(i1)) return
        i2 = this%vert2cell_ja(this%vert2cell_ia(i1) + j - 1)
      end select
    case (IDX_EDGE, IDX_FACE)
      select case (idx2_type)
      case (IDX_VERTEX)
        i2 = this%edge2vert(j,i1)
      case (IDX_EDGE, IDX_FACE)
        if (j > this%edge2edge_ia(i1 + 1) - this%edge2edge_ia(i1)) return
        i2 = this%edge2edge_ja(this%edge2edge_ia(i1) + j - 1)
      case (IDX_CELL)
        i2 = this%edge2cell(j,i1)
        if (i2 == INVALID_NEIGHBOUR) return
      end select
    case (IDX_CELL)
      select case (idx2_type)
      case (IDX_VERTEX)
        i2 = this%cell2vert(j,i1)
      case (IDX_EDGE, IDX_FACE)
        i2 = this%cell2edge(j,i1)
      case (IDX_CELL)
        i2 = this%cell2cell(j,i1)
        if (i2 == INVALID_NEIGHBOUR) return
      end select
    end select

    status = .true.
    idx2(1) = i2
  end subroutine

  subroutine triang_grid_get_adjoint(this, idx, len, surf, vol)
    class(triang_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! cell indices. size: (idx_dim=1)
    real, optional,     intent(out) :: len(:,:)
      !! edge lengths (max_cell_nedge=3, idx_dim=1)
    real, optional,     intent(out) :: surf(:,:)
      !! adjoint surface parts per edge (max_cell_nedge=3, idx_dim=1)
    real, optional,     intent(out) :: vol(:)
      !! adjoint volume parts per vertex (cell_nvert=3)

    integer :: iv1, iv2, ie, ii, jv1, jv2, jv3
    real    :: p(2,3), len_(3), surf_(3), R

    ! get vertex coordinates
    call this%get_cell(idx, p)

    ! get edge lengths
    do ii = 1, 3
      ie  = this%cell2edge(ii,idx(1))
      iv1 = this%edge2vert( 1,ie)
      iv2 = this%edge2vert( 2,ie)

      p(:,1) = this%vert(:,iv1)
      p(:,2) = this%vert(:,iv2)

      len_(ii) = sqrt((p(1,2) - p(1,1))**2 + (p(2,2) - p(2,1))**2)
    end do
    if (present(len)) then
      len(:,1) = len_
    end if

    ! circumscribed radius
    R = (len_(1) * len_(2) * len_(3)) / sqrt((len_(1) + len_(2) + len_(3)) * ( len_(1) + len_(2) - len_(3)) * &
      &                                      (len_(1) - len_(2) + len_(3)) * (-len_(1) + len_(2) + len_(3)))

    ! adjoint surface parts
    do ii = 1, 3
      if (R <= 0.5 * len_(ii)) then
        surf_(ii) = 0
      else
        surf_(ii) = sqrt(R**2 - (0.5 * len_(ii))**2)
      end if
    end do
    if (present(surf)) then
      surf(:,1) = surf_
    end if

    ! adjoint volume parts
    if (present(vol)) then
      vol = 0

      jv1 = this%cell2vert(1,idx(1))
      jv2 = this%cell2vert(2,idx(1))
      jv3 = this%cell2vert(3,idx(1))

      do ii = 1, 3
        ie = this%cell2edge(ii,idx(1))
        iv1 = this%edge2vert(1,ie)
        iv2 = this%edge2vert(2,ie)

        if (iv1 == jv1) then
          vol(1) = vol(1) + 0.5 * len_(ii) * surf_(ii)
        elseif (iv1 == jv2) then
          vol(2) = vol(2) + 0.5 * len_(ii) * surf_(ii)
        else
          vol(3) = vol(3) + 0.5 * len_(ii) * surf_(ii)
        end if

        if (iv2 == jv1) then
          vol(1) = vol(1) + 0.5 * len_(ii) * surf_(ii)
        elseif (iv2 == jv2) then
          vol(2) = vol(2) + 0.5 * len_(ii) * surf_(ii)
        else
          vol(3) = vol(3) + 0.5 * len_(ii) * surf_(ii)
        end if
      end do
    end if
  end subroutine

  subroutine triang_grid_output(this, of, unit)
    !! output triangle grid
    class(triang_grid),     intent(in)    :: this
    type(output_file),      intent(inout) :: of
      !! output file handle
    character(*), optional, intent(in)    :: unit
      !! physical unit of coordinates; default = "um"

    character(:), allocatable  :: unit_
    type(json_object), pointer :: obj

    unit_ = "um"
    if (present(unit)) unit_ = unit

    obj => of%new_object("Grids")
    call obj%add_string("Name", this%name)
    call obj%add_string("Type", "Triangle")
    call of%write(obj, "Vertices", this%vert, unit_)
    call of%write(obj, "Cells", this%cell2vert)
    call of%write(obj, "Edges", this%edge2vert)
  end subroutine

  subroutine triang_grid_output_plotmtv(this, fname, unit)
    !! write grid to plotmtv file
    class(triang_grid), intent(in) :: this
    character(*),       intent(in) :: fname
      !! output file name, e.g. "output/tmp/triang.plt"
    character(*),       intent(in) :: unit
      !! unit of grid variable, e.g. "cm"

    integer               :: i
    type(plotmtv)         :: pmtv
    type(plotset_options) :: opts

    ! init plotmtv file handle
    opts%equalscale = .true.
    call pmtv%init(fname)
    call pmtv%write_header(plotset_opts=opts)

    ! output edges as curves
    associate (e2v => this%edge2vert)
      do i = 1, this%nedge
        call pmtv%write_curve(denorm(this%vert(1,e2v(:,i)), unit), denorm(this%vert(2,e2v(:,i)), unit))
      end do
    end associate

    call pmtv%close()
  end subroutine

  subroutine triang_grid_output_matlab(this, fname, unit)
    !! write grid to plot with matlab's routine triplot(c2v,x,y)
    class(triang_grid), intent(in) :: this
    character(*),       intent(in) :: fname
      !! output file name
    character(*),       intent(in) :: unit
      !! unit of grid variable, e.g. "nm"

    integer :: iounit, iostat, i

    ! vertices
    open (newunit = iounit, file = fname // "_vert.dat", status = "new", action = "write", iostat = iostat)
    if (iostat /= 0) call program_error("File could not be opened")
    do i = 1, this%nvert
      write (iounit, "(ES25.16E3,ES25.16E3)") denorm(this%vert(:,i), unit)
    end do
    close (iounit)

    ! triangle connectivity
    open (newunit = iounit, file = fname // "_c2v.dat", status = "new", action = "write", iostat = iostat)
    if (iostat /= 0) call program_error("File could not be opened")
    do i = 1, this%ncell
      write (iounit, "(I0, x, I0, x, I0)") this%cell2vert(1:3,i)
    end do
    close (iounit)
  end subroutine

  subroutine triang_grid_read_matlab(this, name, fname, unit)
    !! read grid from data files made for plotting with Matlab (see subroutine triang_grid_output_matlab)
    class(triang_grid), intent(out) :: this
    character(*),       intent(in)  :: name
      !! grid name
    character(*),       intent(in)  :: fname
      !! output file name
    character(*),       intent(in)  :: unit
      !! unit of grid variable, e.g. "nm"

    integer              :: iounit, iostat, i, cell2vert(3)
    integer, allocatable :: cells(:,:)
    real                 :: coord(2)
    real,    allocatable :: vert(:,:)
    type(vector_int)     :: c2v(3)
    type(vector_real)    :: v(2)

    ! init vecs
    do i = 1, 3
      if (i < 3) call v(i)%init(0)
      call c2v(i)%init(0)
    end do

    ! read vertices
    open (newunit = iounit, file = fname // "_vert.dat", status = "old", action = "read", iostat = iostat)
    if (iostat /= 0) call program_error("File could not be opened")

    do while (iostat == 0)
      read (iounit, fmt="(ES25.16E3,ES25.16E3)", iostat=iostat) coord
      if (iostat == 0) then
        coord = norm(coord, unit)
        do i = 1, 2
          call v(i)%push(coord(i))
        end do
      end if
    end do
    close (iounit)

    ! read triangle connectivity
    open (newunit = iounit, file = fname // "_c2v.dat", status = "old", action = "read", iostat = iostat)
    if (iostat /= 0) call program_error("File could not be opened")

    do while (iostat == 0)
      read (iounit, fmt=*, iostat=iostat) cell2vert
      if (iostat == 0) then
        do i = 1, 3
          call c2v(i)%push(cell2vert(i))
        end do
      end if
    end do
    close (iounit)

    ! init triangular grid
    allocate (vert(2,v(1)%n), cells(3,c2v(1)%n))
    do i = 1, 3
      if (i < 3) vert(i,:) = v(i)%to_array()
      cells(i,:) = c2v(i)%to_array()
    end do
    call this%init(name, vert, cells)
  end subroutine

end module
