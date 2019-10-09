#!/usr/bin/env python

import os
import sys

for i in range(1, len(sys.argv)):
  name = sys.argv[i]
  filename = "modules/test_"+name+".f90"

  if os.path.exists(filename):
    print("test "+name+" already exists!")
    continue
  else:
    print("generating test "+name)

  f = open(filename, "w+")

  f.write("module test_"+name+"_m\n")
  f.write("  use test_case_m\n")
  f.write("  use "+name+"_m\n")
  f.write("  implicit none\n\n")
  f.write("contains\n\n")
  f.write("  subroutine test_"+name+"()\n")
  f.write("    type(test_case) :: tc\n\n")
  f.write("    print \"(1A)\", \"test_"+name+"\"\n")
  f.write("    call tc%init(\""+name+"\")\n\n\n\n")
  f.write("    call tc%finish()\n")
  f.write("  end subroutine test_"+name+"\n\n")
  f.write("end module test_"+name+"_m\n")

  f.close()

print("generating test.f90")
test_modules = [t.replace("test_","").replace(".f90","") for t in os.listdir("modules/") if os.path.isfile(os.path.join("modules/", t))]

if os.path.exists("test.f90"):
  os.remove("test.f90")

f = open("test.f90", "w+")

f.write("program test\n")
for t in test_modules:
  f.write("  use test_"+t+"_m\n")
f.write("  implicit none\n\n")
f.write("#define all 1\n\n")

for i in range(0, len(test_modules)):
  t = test_modules[i]
  f.write("#define "+t+" "+str(i+2)+"\n")
  f.write("#if ((USE_TEST=="+t+") .or. (USE_TEST==all))\n")
  f.write("  call test_"+t+"()\n")
  f.write("#endif\n")
  f.write("#undef "+t+"\n\n")

f.write("#undef all\n\n")

f.write("end program test\n")