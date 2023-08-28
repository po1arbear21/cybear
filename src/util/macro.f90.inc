m4_divert(-1)

! change m4 quotes to curly brackets
m4_changequote({,})

! suppress unused argument warnings
m4_define({m4_ignore}, {associate( $1 => $1 ); end associate})

! assert macro for error checking, only executed in debug mode
m4_define({m4_assert}, m4_ifdef({m4_debug},{if (.not. ($1)) call assert_failed("$1", "{m4___file__}", {m4___line__})}))

! macros for quick logging
m4_define({m4_info}, call logging(2, {$1 m4_ifelse($#,2,{, }$2)}, file = "{m4___file__}", line = {m4___line__}))
m4_define({m4_error}, call logging(4, {$1 m4_ifelse($#,2,{, }$2)}, file = "{m4___file__}", line = {m4___line__}))

! array shape macro (e.g. m4_shape(4) => :,:,:,:)
m4_define({m4_shape},{m4_ifelse($1,1,:,{m4_shape(m4_decr($1)),:})})
m4_define({m4_pshape},{m4_ifelse($1,0,,{(m4_shape($1))})})

! get expanded typename (e.g. int => integer; my_type => type(my_type))
m4_define({m4_type},{m4_dnl
m4_ifelse($1,int,integer,{m4_dnl
m4_ifelse($1,int32,integer(int32),{m4_dnl
m4_ifelse($1,int64,integer(int64),{m4_dnl
m4_ifelse($1,log,logical,{m4_dnl
m4_ifelse($1,real,real,{m4_dnl
m4_ifelse($1,cmplx,complex,{m4_dnl
m4_ifelse($1,char,character,{m4_dnl
type($1)m4_dnl
})})})})})})})})

! get pointer size measured in integers
m4_define({m4_ptrsize},{m4_ifelse(m4_intsize,64,1,{m4_ifelse(m4_intsize,32,2,0)})})

m4_divert(0)m4_dnl
