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

## Fargo

The project can be built and run using our custom build-system named **fargo**. The configuration is contained in the
file `fargo.toml`. See the **fargo documentation** at `/home/pd/fargo/doc/` for further information.

To compile and run the job specified in the first section in `fargo.toml` you can simply use the command
```bash
$ fargo
```

To run a different section run
```bash
$ fargo run my_section_name
```

To *clean* jobs and build artifacts run
```bash
$ fargo clean
```

If you get strange Python errors, try to delete the `run` and the `build` folders (backup any data you want to keep
first).

## Library Support

The following libraries are currently installed and available at the ITHE:
* ARPACK
* EXPOKIT
* FEAST
* ILUPACK
* MPFR
* MUMPS
* QUADPACK
* SPIKE
* TRIANGLE

The libraries are located at `/home/pd/library`. To be able to use them, simply add the following line to your
`.bashrc`:
```bash
source /home/pd/library/set_lib_paths.sh
```

## Setup on your own computer at home
If you want to use this software on your own computer, you need to clone and compile the **fortran-basic-libs**
repository first. For further instructions see the README.md file of that repository.
