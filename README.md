# Fortran Basic

Fortran project template can be used for your own project.

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

To *clean* build files (with/without library build files, e.g. expokit) run
```bash
$ make clean
$ make clean_all
```

To create *documentation* run
```bash
$ make doc
```

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
