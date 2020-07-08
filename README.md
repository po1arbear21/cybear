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
To compile in (default) *debug* mode run
```bash
$ make
```

To compile in *release* mode (includes optimization) run
```bash
$ make BUILD=release
```

For *fast/parallel* compilation run
```bash
$ make -j4
```

To *clean* build files (with/without library build files) run

```bash
$ make clean
$ make clean_all
```

To create *documentation* run
```bash
$ make doc
```

To use *gfortran* instead of *ifort* run once
```bash
$ make blas95 lapack95 COMPILER=gnu
```
After blas95 and lapack95 libraries have been created run
```bash
$ make COMPILER=gnu
```

### ILUPACK Support
1. If you want to use the ILUPACK solver, first download the zip-archive from http://ilupack.tu-bs.de/
  * Select *Linux 64 BIT version using long integer* and use *MUMPS-based matching*.
  * You can download the *GNU* and/or *INTEL* depending on your compiler choice.
2. Extract the archive and copy the folder *lib/GNU64_long* and/or *lib/INTEL64_long* to a new folder in your home folder name e.g. ILUPACK.
3. Add to your *.bashrc* (replace *YOUR_NAME_HERE* with your name)
  ```bash
  export ILUPACKROOT='/home/YOUR_NAME_HERE/ILUPACK'
  ```

If *ILUPACKROOT* exists, it will automatically be enabled in your project. To disable ILUPACK you can use the option
```bash
$ make USE_ILUPACK=false
```