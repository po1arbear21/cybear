# Fortran Basic

Fortran project template that can be used for your own project.

## Initialization
1. Setup an account at git.rwth-aachen.de.
2. Get invited to the ITHE group.
3. Get invited to this project.
4. Fork this project.
5. Rename fork and change path.
6. Clone the forked project: `$ git clone <url> <folder>`
7. Set up the Fortran Basic repo as a remote such that you can pull new updates from time to time
  * set up remote: `$ git remote add upstream git@git.rwth-aachen.de:ithe/fortran-basic`
  * pull new updates:
```bash
$ git fetch upstream
$ git checkout master
$ git merge upstream/master
```

## Usage
### Makefile
Set default options in *options.mk*.

To compile (with default options) run
```bash
$ make
```

Options can be overwritten by passing values to make, e.g.
```bash
$ make BUILD=release
```

For *fast/parallel* compilation use the -j option:
```bash
$ make -j4
```

To *clean* build files run

```bash
$ make clean
```

To create *documentation* run
```bash
$ make doc
```

### Library Support
The intel MKL must be installed on your system, and ILUPACKROOT should be set as an environment variable (e.g. in ~/.bashrc).
If you want to use any supported additional libraries, clone the git repository fortran-basic-libs and follow the instructions.

Libraries currently contained in *fortran-basic-libs*:
  * ARPACK
  * EXPOKIT
  * FEAST
  * ILUPACK
  * MUMPS
  * QUADPACK
