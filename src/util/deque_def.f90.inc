type deque_{}T
  !! double-ended queue
  m4_type(T), allocatable :: d(:)
    !! data array
  integer :: n
    !! number of elements
  integer :: ifront
    !! front index
  integer :: iback
    !! back index
contains
  procedure :: init       => deque_{}T{}_init
  procedure :: destruct   => deque_{}T{}_destruct
  procedure :: reset      => deque_{}T{}_reset
  procedure :: reserve    => deque_{}T{}_reserve
  procedure :: resize     => deque_{}T{}_resize
  procedure :: front      => deque_{}T{}_front
  procedure :: back       => deque_{}T{}_back
  procedure :: push_front => deque_{}T{}_push_front
  procedure :: push_back  => deque_{}T{}_push_back
  procedure :: pop_front  => deque_{}T{}_pop_front
  procedure :: pop_back   => deque_{}T{}_pop_back
  procedure :: shrink     => deque_{}T{}_shrink
  procedure :: to_array   => deque_{}T{}_to_array
  procedure :: from_array => deque_{}T{}_from_array
end type

m4_undefine({T})
