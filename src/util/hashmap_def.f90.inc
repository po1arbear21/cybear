type hashmap_{}T
  !! hashmap using integer (array) keys
  type(vector_int)     :: keys
    !! key data
  type(vector_int)     :: keys_i0
    !! key start indices
  type(vector_{}T)     :: values
    !! value data
  integer, allocatable :: table(:)
    !! lookup table (used internally)
contains
  procedure :: init     => hashmap_{}T{}_init
  procedure :: destruct => hashmap_{}T{}_destruct
  procedure :: reset    => hashmap_{}T{}_reset
  procedure :: reserve  => hashmap_{}T{}_reserve
  generic   :: get      => hashmap_{}T{}_get_key_1, hashmap_{}T{}_get_key_n
  generic   :: set      => hashmap_{}T{}_set_key_1, hashmap_{}T{}_set_key_n

  procedure, private :: lookup => hashmap_{}T{}_lookup
  procedure, private :: hashmap_{}T{}_get_key_1
  procedure, private :: hashmap_{}T{}_get_key_n
  procedure, private :: hashmap_{}T{}_set_key_1
  procedure, private :: hashmap_{}T{}_set_key_n
end type

m4_undefine({T})
