#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define DATA_SHAPE0
#define DATA_SHAPE1 :
#define DATA_SHAPE2 :,:
#define DATA_SHAPE3 :,:,:
#define DATA_SHAPE4 :,:,:,:
#define DATA_SHAPE5 :,:,:,:,:
#define DATA_SHAPE6 :,:,:,:,:,:
#define DATA_SHAPE7 :,:,:,:,:,:,:
#define DATA_SHAPE8 :,:,:,:,:,:,:,:
#define DATA_SHAPE(N) DATA_SHAPE ## N

#define OUTPUT_FILE_WRITE    CONCAT3(output_file_write, T, N)


integer :: asdf(DATA_SHAPE(4))

#undef T
#undef TT
#undef N

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef DATA_SHAPE0
#undef DATA_SHAPE1
#undef DATA_SHAPE2
#undef DATA_SHAPE3
#undef DATA_SHAPE4
#undef DATA_SHAPE5
#undef DATA_SHAPE6
#undef DATA_SHAPE7
#undef DATA_SHAPE8
#undef DATA_SHAPE

#undef OUTPUT_FILE_WRITE

! subroutine OUTPUT_FILE_WRITE(this, obj, name, values)
!   !! write data to binary file
!   class(output_file),         intent(inout) :: this
!   type(json_object), pointer, intent(inout) :: obj
!     !! parent json object
!   character(*),               intent(in)    :: name
!     !! data name
!   TT,                         intent(in)    :: values(DATA_SHAPE)
!     !! data array
! end subroutine


! subroutine output_file_write_i32_2(this, obj, name, values)
!   !! write 2D 32-bit integer array to binary file
!   class(output_file),         intent(inout) :: this
!   type(json_object), pointer, intent(inout) :: obj
!     !! parent json object
!   character(*),               intent(in)    :: name
!     !! data name
!   integer(kind=int32),        intent(in)    :: values(:,:)
!     !! 2D 32-bit integer array

!   integer                    :: i
!   integer, allocatable       :: tmp(:)
!   type(json_array),  pointer :: sh
!   type(json_object), pointer :: dat

!   ! create json shape array
!   allocate (sh)
!   call sh%init()
!   tmp = shape(values)
!   do i = 1, size(tmp)
!     call sh%add(tmp(i))
!   end do

!   ! create json data descriptor and add it to parent object
!   allocate (dat)
!   call dat%init()
!   call dat%add("Type", "i32")
!   call dat%add("Shape", sh)
!   call dat%add("Index", this%index_bin)
!   call obj%add(name, dat)

!   ! write data to binary file
!   write (this%funit_bin) values
!   this%index_bin = this%index_bin + size(values) * 4
! end subroutine
