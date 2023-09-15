from typing import List, Literal
import os

class Config:
    FC:           str       = "ifort"
    FFLAGS:       List[str] = ["-warn", "all", "-qopenmp", "-fp-model", "precise", "-real-size", "64", "-fpp"] # No '-march=native'!
    LFLAGS:       List[str] = []#"-Wl,--no-warn-execstack"]
    FMODULE:      str       = "-module"
    FSYNTAX_ONLY: str       = "-syntax-only"

    CC:     str       = "icc"
    CFLAGS: List[str] = []

    LIBS: List[str] = []
    FINC: List[str] = []
    CINC: List[str] = []

    def __init__(self, mode: Literal["debug", "profile", "release"], intsize: Literal[32, 64], idxsize: Literal[32, 64], libraries: List[str]) -> None:

        if mode == "debug":
            self.FFLAGS += ["-O0", "-qmkl=sequential", "-g", "-check", "all", "-check", "noarg_temp_created", "-fpe1", "-traceback", "-debug", "extended", "-init=huge"]
            self.CFLAGS += ["-O0", "-g"]
        elif mode == "profile":
            self.FFLAGS += ["-O3", "-qmkl", "-ftz", "-g", "-shared-intel", "-debug", "inline-debug-info", "-parallel-source-info=2"]
            self.CFLAGS += ["-O3", "-g", "-shared-intel"]
        elif mode == "release":
            self.FFLAGS += ["-O3", "-qmkl", "-ftz"]
            self.CFLAGS += ["-O3"]

        if intsize == 64:
            self.FFLAGS += ["-i8"]

        def get_env(key: str) -> str:
            if value := os.getenvb(key.encode(), b"").decode():
                return value
            quit()

        def addlib(*path_parts) -> None:
            self.LIBS.append(os.path.join(*path_parts))

        arch = "lp64" if intsize == 32 else "ilp64"

        # FEAST
        if "FEAST" in libraries:
            feastroot = get_env("FEASTROOT")
            addlib(feastroot, "lib", f"intel_{arch}", "libfeast.a")

        # SPIKE
        if "SPIKE" in libraries:
            spikeroot = get_env("SPIKEROOT")
            addlib(spikeroot, "lib", f"intel_{arch}", "libspike.a")

        # BLAS95 and LAPACK95
        mklroot = get_env("MKLROOT")

        addlib(mklroot, "lib", "intel64", f"libmkl_blas95_{arch}.a")
        addlib(mklroot, "lib", "intel64", f"libmkl_lapack95_{arch}.a")

        # ARPACK
        if "ARPACK" in libraries:
            arpackroot = get_env("ARPACKROOT")
            self.LIBS.append(
                os.path.join(arpackroot, "lib", f"intel_{arch}", "libarpack.a")
            )

        # EXPOKIT
        if "EXPOKIT" in libraries:
            expokitroot = get_env("EXPOKITROOT")
            addlib(expokitroot, "lib", f"intel_{arch}", "expokit.a")

        # ILUPACK
        if "ILUPACK" in libraries:
            ilupackroot = get_env("ILUPACKROOT")
            archilu = "ilp64" if idxsize == 64 else "lp64"

            self.LIBS.append("-Wl,--start-group")
            addlib(ilupackroot, "lib", f"intel_{archilu}", "libilupack_mumps.a")
            addlib(ilupackroot, "lib", f"intel_{archilu}", "libamd.a")
            addlib(ilupackroot, "lib", f"intel_{archilu}", "libblaslike.a")
            addlib(ilupackroot, "lib", f"intel_{archilu}", "libcamd.a")
            addlib(ilupackroot, "lib", f"intel_{archilu}", "libmetis.a")
            addlib(ilupackroot, "lib", f"intel_{archilu}", "libmetisomp.a")
            addlib(ilupackroot, "lib", f"intel_{archilu}", "libsparspak.a")
            addlib(ilupackroot, "lib", f"intel_{archilu}", "libsuitesparseconfig.a")
            addlib(ilupackroot, "lib", f"intel_{archilu}", "libmumps.a")
            self.LIBS.append("-Wl,--end-group")

        # MUMPS
        if "MUMPS" in libraries:
            mumpsroot = get_env("MUMPSROOT")
            metisroot = get_env("METISROOT")
            scotchroot = get_env("SCOTCHROOT")

            self.LIBS.append("-Wl,--start-group")
            addlib(mumpsroot, "lib", f"intel_{arch}", "libdmumps.a")
            addlib(mumpsroot, "lib", f"intel_{arch}", "libzmumps.a")
            addlib(mumpsroot, "lib", f"intel_{arch}", "libmumps_common.a")
            addlib(mumpsroot, "lib", f"intel_{arch}", "libpord.a")
            addlib(mumpsroot, "lib", f"intel_{arch}", "libmpiseq.a")
            addlib(metisroot, "lib", "intel_ilp64", "libmetis.a")
            addlib(scotchroot, "lib", "intel_ilp64", "libesmumps.a")
            addlib(scotchroot, "lib", "intel_ilp64", "libscotch.a")
            addlib(scotchroot, "lib", "intel_ilp64", "libscotcherr.a")
            self.LIBS.append("-Wl,--end-group")

            self.FINC.append("-I" + os.path.join(mumpsroot, "include", f"intel_{arch}"))

        # QUADPACK
        if "QUADPACK" in libraries:
            quadpackroot = get_env("QUADPACKROOT")
            addlib(quadpackroot, "lib", f"intel_{arch}", "quadpack.a")

        # MPFR
        if "MPFR" in libraries:
            gmproot = get_env("GMPROOT")
            mpfrroot = get_env("MPFRROOT")

            addlib(mpfrroot, "lib", "intel_ilp64", "libmpfr.a")
            addlib(gmproot, "lib", "intel_ilp64", "libgmp.a")
            self.CINC.append("-I" + os.path.join(mpfrroot, "include", "intel_ilp64"))
            self.CINC.append("-I" + os.path.join(gmproot, "include", "intel_ilp64"))

        # TRIANGLE
        if "TRIANGLE" in libraries:
            triangleroot = get_env("TRIANGLEROOT")
            addlib(triangleroot, "lib", "intel_lp64", "triangle.a")
