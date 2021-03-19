module test_hashmap_m
  use array_m
  use hashmap_m
  use test_case_m

  implicit none

  private
  public test_hashmap

contains

  subroutine test_hashmap()
    type(test_case) :: tc

    call tc%init("hashmap")

    block
      integer           :: i, value
      logical           :: status
      type(hashmap_int) :: hmap

      ! init
      call hmap%init(c=32)
      call tc%assert(allocated(hmap%keys%d  ), "hashmap_int init: allocated(hmap%keys%d)")
      call tc%assert(allocated(hmap%values%d), "hashmap_int init: allocated(hmap%values%d)")
      call tc%assert(allocated(hmap%table   ), "hashmap_int init: allocated(hmap%table)")
      call tc%assert_eq(32, size(hmap%keys%d)  , "hashmap_int init: keys capacity")
      call tc%assert_eq(32, size(hmap%values%d), "hashmap_int init: values capacity")

      ! set/get
      call hmap%init(c = 32)
      call hmap%set(32354, 411)
      call hmap%set(-47882, 923)
      call hmap%set(12234, 158)
      call hmap%set(12234, 159)
      call hmap%get(32354, value)
      call tc%assert_eq(411, value, "hashmap_int set/get 1")
      call hmap%get(-47882, value)
      call tc%assert_eq(923, value, "hashmap_int set/get 2")
      call hmap%get(12234, value)
      call tc%assert_eq(159, value, "hashmap_int set/get 3")
      call hmap%get(12235, value, status)
      call tc%assert(.not. status, "hashmap_int set/get 4")

      ! set/get many elements
      do i = 1, 1024
        call hmap%set(i, 2*i)
      end do
      do i = 1, 1024
        call hmap%get(i, value)
        call tc%assert_eq(2*i, value, "hashmap_int set/get many")
      end do

      ! reset
      call hmap%reset()
      call tc%assert_eq(0, hmap%keys%n, "hashmap_int reset: keys%n == 0")
      call tc%assert_eq(0, hmap%values%n, "hashmap_int reset: values%n == 0")
      call tc%assert(all(hmap%table == -1), "hashmap_int reset: all(table == -1)")

      ! destruct
      call hmap%destruct()
      call tc%assert(.not. allocated(hmap%keys%d  ), "hashmap_int destruct: deallocated(hmap%keys%d)")
      call tc%assert(.not. allocated(hmap%values%d), "hashmap_int destruct: deallocated(hmap%values%d)")
      call tc%assert(.not. allocated(hmap%table   ), "hashmap_int destruct: deallocated(hmap%table)")
    end block

    block
      type(hashmap_int) :: hmap
      integer           :: key(2), value
      logical           :: status

      ! set/get
      call hmap%init(keysize=2, c=32)
      key = [4124, 142111]
      call hmap%set(key, 12345)
      key = [5629, 666]
      call hmap%set(key, 4444)
      key = [4124, 143111]
      call hmap%set(key, 54321)
      key = [5629, 666]
      call hmap%set(key, 3333)

      key = [4124, 142111]
      call hmap%get(key, value)
      call tc%assert_eq(12345, value, "hashmap2_int set/get 1")
      key = [4124, 143111]
      call hmap%get(key, value)
      call tc%assert_eq(54321, value, "hashmap2_int set/get 2")
      key = [5629, 666]
      call hmap%get(key, value)
      call tc%assert_eq(3333, value, "hashmap2_int set/get 3")
      key = [123, 456]
      call hmap%get(key, value, status)
      call tc%assert(.not. status, "hashmap2_int set/get 4")
    end block

    call tc%finish()
  end subroutine

end module
