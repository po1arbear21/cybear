type vector_{}T
  m4_type(T), allocatable :: d(:)
    !! data array
  integer         :: n = 0
    !! size of vector
contains
  procedure :: init       => vector_{}T{}_init
  procedure :: destruct   => vector_{}T{}_destruct
  procedure :: reset      => vector_{}T{}_reset
  procedure :: reserve    => vector_{}T{}_reserve
  procedure :: resize     => vector_{}T{}_resize
  procedure :: front      => vector_{}T{}_front
  procedure :: back       => vector_{}T{}_back
  generic   :: push       => vector_{}T{}_push, vector_{}T{}_push_n
  procedure :: pop        => vector_{}T{}_pop
  procedure :: shrink     => vector_{}T{}_shrink
  procedure :: to_array   => vector_{}T{}_to_array
  procedure :: from_array => vector_{}T{}_from_array

  procedure, private :: vector_{}T{}_push
  procedure, private :: vector_{}T{}_push_n
end type

m4_undefine({T})
