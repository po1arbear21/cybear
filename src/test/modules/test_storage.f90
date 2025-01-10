m4_include(../../util/macro.f90.inc)

module test_storage_m

  use test_case_m, only: test_case
  use storage_m,   only: storage, STORAGE_READ, STORAGE_WRITE, DYNAMIC_APP, DYNAMIC_EXT, COMPR_ZLIB, COMPR_BLOSC
  use string_m
  
  implicit none

  private
  public test_storage

contains

  subroutine test_storage()

    type(test_case) :: tc
    type(storage) :: store, shop
    real :: varr, rarr(10), dynrr(5,3)
    integer(kind=8) :: vari,iarr(3,60,200)
    complex :: varc
    logical :: varl
    type(string) :: str
    type(string), allocatable :: str_arr(:)
    real :: dynamic_append(15)

    integer :: dyn(3)
    integer, allocatable :: extract(:)

    real, allocatable :: trarr(:), dyneee(:,:), dynamic_append_return(:)
    integer(kind=8), allocatable :: tiarr(:,:,:)
    
    call tc%init("storage")

    call store%open("test_storage.fbs", STORAGE_WRITE)

    dyn = [1, 2, 3]
    call store%write("dynamic", dyn(1), dynamic=DYNAMIC_EXT)
    call store%write("dynamic", dyn(2), dynamic=DYNAMIC_EXT)
    call store%write("dynamic", dyn(3), dynamic=DYNAMIC_EXT)

    iarr = 2
    m4_ifdef({m4_zlib},{call store%write(new_string("iarr_compressed"), iarr, compression=COMPR_ZLIB)})
    m4_ifdef({m4_blosc},{call store%write(new_string("iarr_compressed_blosc"), iarr, compression=COMPR_BLOSC)})
    call store%write(new_string("iarr"), iarr)

    
    dynamic_append = [-0.8, 7.1, 4.9732, 45e12, -25e-32, 90.3, 1.0, 0.0, 123456.789, 1.9, 1.0, 1.111111111111, 1.0999999999999, -0.0000000001, 0.0000000001]
    call store%write("dynamic_append", dynamic_append(1:6), dynamic=DYNAMIC_APP)

    dynrr(:,1) = [1.0, 2.0, 3.0, 4.0, 5.0]
    call store%write("dynr", dynrr(:,1), dynamic=DYNAMIC_EXT)

    varr = 3.141529
    call store%write(new_string("var"), varr)
    varl = .true.
    call store%write(new_string("logical"), varl)
    varr = 1.602e-19
    call store%write(new_string("q"), varr)
    varc = complex(2.0, -3.96719047e12)
    call store%write(new_string("comp"), varc)

    call store%write(new_string("abccccccc"), [new_string("a"), new_string("bb"), new_string("ccc"), new_string("dddd"), new_string("eeeee"), new_string("ffffff")])

    dynrr(:,2) = [1.4, 2.4, 3.4, 4.4, 5.4]
    call store%write(new_string("dynr"), dynrr(:,2), dynamic=DYNAMIC_EXT)

    call store%write("dynamic_append", dynamic_append(7:15), dynamic=DYNAMIC_APP)

    vari = 3482344
    call store%write(new_string("inttt"), vari)

    dynrr(:,3) = [1.1, 2.1, 3.1, 4.1, 5.1]
    call store%write(new_string("dynr"), dynrr(:,3), dynamic=DYNAMIC_EXT)

    rarr = [3.290, 0.433, 0.77234e-12, 0.87653, 0.72, 47.5e33, 0.99999999, 1.0, 0.9e99, 1.9]
    call store%write(new_string("rarr"), rarr)

    call store%write(new_string("string"), new_string("This is a test string for ssssssssssssss."))

    call store%close()
    call shop%open("test_storage.fbs", STORAGE_READ)

    call shop%read(new_string("dynamic"), extract)
    call tc%assert(allocated(extract), "Did not allocate array")
    call tc%assert_eq(dyn, extract, "Retrieved dynamic array does not match input")

    call shop%read(new_string("inttt"), vari)
    call tc%assert_eq(3482344, int(vari), "Retrieved variable inttt does not match input")
    call shop%read(new_string("var"), varr)
    call tc%assert_eq(3.141529, varr, 1e-16, 1e-16, "Retrieved variable var does not match input")
    call shop%read(new_string("logical"), varl)
    call tc%assert(varl, "Retrieved variable logical does not match input")
    call shop%read(new_string("q"), varr)
    call tc%assert_eq(1.602e-19, varr, 1e-16, 1e-16, "Retrieved variable q does not match input")
    call shop%read(new_string("comp"), varc)
    call tc%assert_eq(complex(2.0, -3.96719047e12), varc, 1e-16, 1e-16, "Retrieved variable comp does not match input")

    call shop%read(new_string("rarr"), trarr)
    call tc%assert(allocated(trarr), "Did not allocate array")
    call tc%assert_eq(rarr, trarr, 1e-16, 1e-16, "Retrieved array trarr does not match")
    call shop%read(new_string("iarr"), tiarr)
    call tc%assert(allocated(tiarr), "Did not allocate array")
    call tc%assert_eq(int(iarr), int(tiarr), "Retrieved array tiarr does not match")
    
    m4_ifdef({m4_zlib},{deallocate(tiarr)
    call shop%read(new_string("iarr_compressed"), tiarr)
    call tc%assert(allocated(tiarr), "Did not allocate array")
    call tc%assert_eq(int(iarr), int(tiarr), "Retrieved array tiarr does not match")})

    m4_ifdef({m4_blosc},{deallocate(tiarr)
    call shop%read(new_string("iarr_compressed_blosc"), tiarr)
    call tc%assert(allocated(tiarr), "Did not allocate array")
    call tc%assert_eq(int(iarr), int(tiarr), "Retrieved array tiarr does not match")})

    call shop%read("dynamic_append", dynamic_append_return)
    call tc%assert_eq(dynamic_append, dynamic_append_return, 1e-16, 1e-16, "The retrieved dynamic append array does match its input")

    call shop%read(new_string("string"), str)
    call tc%assert(new_string("This is a test string for ssssssssssssss.") == str, "Strings do not match")

    call shop%read(new_string("dynr"), dyneee)
    call tc%assert(allocated(trarr), "Did not allocate array")
    call tc%assert_eq(dynrr, dyneee, 1e-16, 1e-16, "Retrieved array does not match")

    call shop%read(new_string("abccccccc"), str_arr)
    call tc%assert_eq([new_string("a"), new_string("bb"), new_string("ccc"), new_string("dddd"), new_string("eeeee"), new_string("ffffff")], str_arr, "String arrays do not match")

    call shop%close()
  
    call tc%finish()
  end

end module