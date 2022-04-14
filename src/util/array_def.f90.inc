m4_ifelse(N,0,{
  type ptr_{}T
    m4_type(T), pointer :: p => null()
  end type
},{
  type array{}N{}_{}T
    m4_type(T), allocatable :: d(m4_shape(N))
  end type

  type ptr{}N{}_{}T
    m4_type(T), pointer :: p(m4_shape(N)) => null()
  end type
})

m4_undefine({N},{T})
