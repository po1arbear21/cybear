module test_qsort_m
  use test_case_m
  use qsort_m
  implicit none

contains

  subroutine test_qsort()
    type(test_case) :: tc

    print "(A)", "test_qsort"
    call tc%init("qsort")

    ! base case: 2 elements
    block
      integer :: ar(2), perm(2)

      ar = [1, 9]
      call qsort(ar, perm = perm)
      call tc%assert_eq([1, 9], ar  , "sort 2 elements (ar); test 1")
      call tc%assert_eq([1, 2], perm, "sort 2 elements (perm); test 1")

      ar = [9, 1]
      call qsort(ar, perm = perm)
      call tc%assert_eq([1, 9], ar  , "sort 2 elements (ar); test 2")
      call tc%assert_eq([2, 1], perm, "sort 2 elements (perm); test 2")
    end block

    ! base case: 3 elements
    block
      integer :: ar(3), perm(3)

      ar = [1, 7, 9]
      call qsort(ar, perm = perm)
      call tc%assert_eq([1, 7, 9], ar  , "sort 3 elements (ar); test 1")
      call tc%assert_eq([1, 2, 3], perm, "sort 3 elements (perm); test 1")

      ar = [9, 1, 7]
      call qsort(ar, perm = perm)
      call tc%assert_eq([1, 7, 9], ar  , "sort 3 elements (ar); test 2")
      call tc%assert_eq([2, 3, 1], perm, "sort 3 elements (perm); test 2")
    end block

    ! base case: 4 elements
    block
      integer :: ar(4), perm(4)

      ar = [1, 7, 8, 9]
      call qsort(ar, perm = perm)
      call tc%assert_eq([1, 7, 8, 9], ar  , "sort 4 elements (ar); test 1")
      call tc%assert_eq([1, 2, 3, 4], perm, "sort 4 elements (perm); test 1")

      ar = [9, 1, 7, 8]
      call qsort(ar, perm = perm)
      call tc%assert_eq([1, 7, 8, 9], ar  , "sort 4 elements (ar); test 2")
      call tc%assert_eq([2, 3, 4, 1], perm, "sort 4 elements (perm); test 2")

      ar = [7, 1, 9, 8]
      call qsort(ar, perm = perm)
      call tc%assert_eq([1, 7, 8, 9], ar  , "sort 4 elements (ar); test 3")
      call tc%assert_eq([2, 1, 4, 3], perm, "sort 4 elements (perm); test 3")
    end block

    ! base case: 5 elements
    block
      integer :: ar(5), perm(5)

      ar = [1, 5, 7, 8, 9]
      call qsort(ar, perm = perm)
      call tc%assert_eq([1, 5, 7, 8, 9], ar  , "sort 5 elements (ar); test 1")
      call tc%assert_eq([1, 2, 3, 4, 5], perm, "sort 5 elements (perm); test 1")

      ar = [8, 9, 1, 5, 7]
      call qsort(ar, perm = perm)
      call tc%assert_eq([1, 5, 7, 8, 9], ar  , "sort 5 elements (ar); test 2")
      call tc%assert_eq([3, 4, 5, 1, 2], perm, "sort 5 elements (perm); test 2")

      ar = [5, 1, 9, 7, 8]
      call qsort(ar, perm = perm)
      call tc%assert_eq([1, 5, 7, 8, 9], ar  , "sort 5 elements (ar); test 3")
      call tc%assert_eq([2, 1, 4, 5, 3], perm, "sort 5 elements (perm); test 3")
    end block

    ! 32 elements
    block
      integer :: ar(32), perm(32), e(32)

      ! no double elements
      ar = [54,89,38,37,100,82,11,17,5,56,84,6,15,26,69,75,34,22,44,94,98,46,71,93,55,90,64,48,28,80,87,57]
      call qsort(ar, perm = perm)
      e = [5,6,11,15,17,22,26,28,34,37,38,44,46,48,54,55,56,57,64,69,71,75,80,82,84,87,89,90,93,94,98,100]
      call tc%assert_eq(e, ar, "sort 32 elements (ar); test 1")
      e = [9,12,7,13,8,18,14,29,17,4,3,19,22,28,1,25,10,32,27,15,23,16,30,6,11,31,2,26,24,20,21,5]
      call tc%assert_eq(e, perm, "sort 32 elements (perm); test 1")

      ar = [50,54,48,89,46,71,23,56,36,29,92,43,19,78,45,2,39,100,52,35,31,44,27,91,90,37,28,60,74,95,14,47]
      call qsort(ar, perm = perm)
      e = [2,14,19,23,27,28,29,31,35,36,37,39,43,44,45,46,47,48,50,52,54,56,60,71,74,78,89,90,91,92,95,100]
      call tc%assert_eq(e, ar, "sort 32 elements (ar); test 2")
      e = [16,31,13,7,23,27,10,21,20,9,26,17,12,22,15,5,32,3,1,19,2,8,28,6,29,14,4,25,24,11,30,18]
      call tc%assert_eq(e, perm, "sort 32 elements (perm); test 2")

      ar = [50,54,50,31,46,54,23,56,36,29,46,43,19,78,45,78,54,50,52,35,31,44,27,91,31,37,31,60,74,95,14,31]
      call qsort(ar)
      e = [14,19,23,27,29,31,31,31,31,31,35,36,37,43,44,45,46,46,50,50,50,52,54,54,54,56,60,74,78,78,91,95]
      call tc%assert_eq(e, ar, "sort 32 elements with duplicates; test 3")
    end block

    ! 37 real elements
    block
      real :: ar(37), e(37)

      ar = [0.8147, 0.9058, 0.1270, 0.9134, 0.6324, 0.0975, 0.2785, 0.5469, &
            0.9575, 0.9649, 0.1576, 0.9706, 0.9572, 0.4854, 0.8003, 0.1419, &
            0.4218, 0.9157, 0.7922, 0.9595, 0.6557, 0.0357, 0.8491, 0.9340, &
            0.6787, 0.7577, 0.7431, 0.3922, 0.6555, 0.1712, 0.7060, 0.0318, &
            0.2769, 0.0462, 0.0971, 0.8235, 0.6948]
      call qsort(ar)
      e = [ 0.0318, 0.0357, 0.0462, 0.0971, 0.0975, 0.1270, 0.1419, 0.1576, &
            0.1712, 0.2769, 0.2785, 0.3922, 0.4218, 0.4854, 0.5469, 0.6324, &
            0.6555, 0.6557, 0.6787, 0.6948, 0.7060, 0.7431, 0.7577, 0.7922, &
            0.8003, 0.8147, 0.8235, 0.8491, 0.9058, 0.9134, 0.9157, 0.9340, &
            0.9572, 0.9575, 0.9595, 0.9649, 0.9706]
      call tc%assert_eq(e, ar, 1e-16, "sort 37 real elements")
    end block

    ! 17 strings
    block
      type(string) :: s(17), e(17)

      s = [ string("nsthemzb"), string("nezlbflp"), string("shjfsiix"), string("jsmyblzt"), &
            string("adghtypa"), string("tquyqeyp"), string("okepmhsw"), string("oxuplhxq"), &
            string("fgityclt"), string("rowvqngo"), string("fqnfgoau"), string("tiewijrx"), &
            string("oyfvgytz"), string("ymutpqwx"), string("gwisvfyc"), string("odhwtxgg"), &
            string("ugpwhttk") ]
      call qsort(s)
      e = [ string("adghtypa"), string("fgityclt"), string("fqnfgoau"), string("gwisvfyc"), &
            string("jsmyblzt"), string("nezlbflp"), string("nsthemzb"), string("odhwtxgg"), &
            string("okepmhsw"), string("oxuplhxq"), string("oyfvgytz"), string("rowvqngo"), &
            string("shjfsiix"), string("tiewijrx"), string("tquyqeyp"), string("ugpwhttk"), &
            string("ymutpqwx") ]
      call tc%assert_eq(e, s, "sort 17 strings")
    end block

    call tc%finish()
  end subroutine

end module
