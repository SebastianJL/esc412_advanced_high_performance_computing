# esc412_advanced_high_performance_computing

Code written by me can be found in src/.


## Build the project

First build blitz.
```bash
cmake -DCMAKE_INSTALL_PREFIX=./ -B ./blitz-build ./blitz-1.0.2
make -C ./blitz-build
make -C ./blitz-build install
```
This will create the folders `include`, `lib` and `blitz-build`. They need to be present in order for the
project to run.

The project can be built with one of the targets in `Makefile` or `Makefile.daint`. The binaries are built into
`make-build`.
Slurm scripts to run the compiled code on daint (or eiger) can be found in scripts/. They are named
`run_<name_of_executable>`.

The cmake stuff can largely be ignored as it does not work and is only present such that the clion code introspection
features work.