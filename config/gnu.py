from typing import List, Literal
import os, sys

class Config:
    FC:           str       = "gfortran"
    FFLAGS:       List[str] = ["-ffree-line-length-none", "-march=native", "-Wall", "-Wno-maybe-uninitialized", "-fopenmp", "-fuse-ld=bfd", "-fdefault-real-8", "-cpp"]
    FMODULE:      str       = "-J"
    FSYNTAX_ONLY: str       = "-fsyntax-only" 

    CC:     str       = "gcc"
    CFLAGS: List[str] = ["-march=native"]

    LIBS: List[str] = []
    FINC: List[str] = []
    CINC: List[str] = []

    def __init__(self, mode: Literal["debug", "profile", "release"], intsize: Literal[32, 64], idxsize: Literal[32, 64], libraries: List[str]) -> None:

        if mode == "debug":
            self.FFLAGS += ["-O0", "-g3", "-ggdb", "-fcheck=all", "-fbacktrace", "-finit-integer=-9999", "-finit-real=-inf"]
            self.CFLAGS += ["-O0", "-g"]
        elif mode == "profile":
            self.FFLAGS += ["-O2", "-g", "-shared-libgcc", "-flto"]
            self.CFLAGS += ["-O3", "-g", "-shared-libgcc", "-flto"]
        elif mode == "release":
            self.FFLAGS += ["-O2", "-flto"]
            self.CFLAGS += ["-O3", "-flto"]

        if intsize == 64:
            self.FFLAGS += ["-fdefault-integer-8"]

        def get_env(key: str) -> str:
            if value := os.getenvb(key.encode(), b"").decode():
                return value
            sys.exit(1)

        def addlib(*path_parts) -> None:
            self.LIBS.append(os.path.join(*path_parts))

        arch = "lp64" if intsize == 32 else "ilp64"

        # FEAST
        if "FEAST" in libraries:
            feastroot = get_env("FEASTROOT")
            addlib(feastroot, "lib", f"gnu_{arch}", "libfeast.a")

        # BLAS95 and LAPACK95
        blasroot = get_env("BLAS95ROOT")
        lapackroot = get_env("LAPACK95ROOT")

        addlib(blasroot, "lib", f"gnu_{arch}", f"libmkl_blas95_{arch}.a")
        addlib(lapackroot, "lib", f"gnu_{arch}", f"libmkl_lapack95_{arch}.a")
        self.FINC.append("-I" + os.path.join(blasroot, "include", f"gnu_{arch}"))
        self.FINC.append("-I" + os.path.join(lapackroot, "include", f"gnu_{arch}"))

        mklroot = get_env("MKLROOT")
        self.FINC.append("-I" + os.path.join(mklroot, "include"))

        self.LIBS.append("-Wl,--start-group")
        addlib(mklroot, "lib", "intel64", f"libmkl_gf_{arch}.a")

        if mode == "debug":
            addlib(mklroot, "lib", "intel64", "libmkl_sequential.a")
        else:
            addlib(mklroot, "lib", "intel64", "libmkl_gnu_thread.a")
        addlib(mklroot, "lib", "intel64", "libmkl_core.a")
        self.LIBS.append("-Wl,--end-group")

        # ARPACK
        if "ARPACK" in libraries:
            arpackroot = get_env("ARPACKROOT")
            self.LIBS.append(
                os.path.join(arpackroot, "lib", f"gnu_{arch}", "libarpack.a")
            )

        # EXPOKIT
        if "EXPOKIT" in libraries:
            expokitroot = get_env("EXPOKITROOT")
            addlib(expokitroot, "lib", f"gnu_{arch}", "expokit.a")

        # ILUPACK
        if "ILUPACK" in libraries:
            ilupackroot = get_env("ILUPACKROOT")
            archilu = "ilp64" if idxsize == 64 else "lp64"

            self.LIBS.append("-Wl,--start-group")
            addlib(ilupackroot, "lib", f"gnu_{archilu}", "libilupack_mumps.a")
            addlib(ilupackroot, "lib", f"gnu_{archilu}", "libamd.a")
            addlib(ilupackroot, "lib", f"gnu_{archilu}", "libblaslike.a")
            addlib(ilupackroot, "lib", f"gnu_{archilu}", "libcamd.a")
            addlib(ilupackroot, "lib", f"gnu_{archilu}", "libmetis.a")
            addlib(ilupackroot, "lib", f"gnu_{archilu}", "libmetisomp.a")
            addlib(ilupackroot, "lib", f"gnu_{archilu}", "libsparspak.a")
            addlib(ilupackroot, "lib", f"gnu_{archilu}", "libsuitesparseconfig.a")
            addlib(ilupackroot, "lib", f"gnu_{archilu}", "libmumps.a")
            self.LIBS.append("-Wl,--end-group")

        # MUMPS
        if "MUMPS" in libraries:
            mumpsroot = get_env("MUMPSROOT")
            metisroot = get_env("METISROOT")
            scotchroot = get_env("SCOTCHROOT")

            self.LIBS.append("-Wl,--start-group")
            addlib(mumpsroot, "lib", f"gnu_{arch}", "libdmumps.a")
            addlib(mumpsroot, "lib", f"gnu_{arch}", "libzmumps.a")
            addlib(mumpsroot, "lib", f"gnu_{arch}", "libmumps_common.a")
            addlib(mumpsroot, "lib", f"gnu_{arch}", "libpord.a")
            addlib(mumpsroot, "lib", f"gnu_{arch}", "libmpiseq.a")
            addlib(metisroot, "lib", "gnu_ilp64", "libmetis.a")
            addlib(scotchroot, "lib", "gnu_ilp64", "libesmumps.a")
            addlib(scotchroot, "lib", "gnu_ilp64", "libscotch.a")
            addlib(scotchroot, "lib", "gnu_ilp64", "libscotcherr.a")
            self.LIBS.append("-Wl,--end-group")

            mkllibs = get_env("MKLLIBS")
            self.LIBS.append(mkllibs)

            self.FINC.append("-I" + os.path.join(mumpsroot, "include", f"gnu_{arch}"))

        # QUADPACK
        if "QUADPACK" in libraries:
            quadpackroot = get_env("QUADPACKROOT")
            addlib(quadpackroot, "lib", f"gnu_{arch}", "quadpack.a")

        # MPFR
        if "MPFR" in libraries:
            gmproot = get_env("GMPROOT")
            mpfrroot = get_env("MPFRROOT")

            addlib(mpfrroot, "lib", "gnu_ilp64", "libmpfr.a")
            addlib(gmproot, "lib", "gnu_ilp64", "libgmp.a")
            self._cinc.append("-I" + os.path.join(mpfrroot, "include", "gnu_ilp64"))
            self._cinc.append("-I" + os.path.join(gmproot, "include", "gnu_ilp64"))

        # TRIANGLE
        if "TRIANGLE" in libraries:
            triangleroot = get_env("TRIANGLEROOT")
            addlib(triangleroot, "lib", "gnu_lp64", "triangle.a")

        self.LIBS.append("-lgomp")
        self.LIBS.append("-lpthread")
        self.LIBS.append("-lm")
        self.LIBS.append("-ldl")