# TDT4200 Problem set 3: Pthreads and OpenMP

## Finite difference approximation of the 2D shallow water equations using Pthreads and OpenMP
In this assignment you will work on an implementation of the Finite Difference Method (FDM) for solving the 2D shallow water equations using POSIX Threads (pthreads) and OpenMP. Information on solving the shallow water equations using FDM and on parallelization with pthreads/OpenMP is described in the lecture slides.

The serial solution can be found in `shallow_water_serial.c` and should be kept as a reference. A skeleton for your parallel implementations can be found in `shallow_water_parallel_{omp, pthreads}.c`, in which you should write your parallel implementations. You should complete the parallel implementation as described by the problem set description.

## Run
### Setup
`make setup`

Creates folders `data`, `plots` and `video`.
- `data`: contains output from the simulation
- `plots`: contains output from plotting
- `video`: contains output from video generation

Compiles the code for comparing solutions.

### Serial solution
**Compile**

`make serial`

**Run**

`./serial -n [grid_size] -i [max_iteration] -s [snapshot_frequency]`

**Example**

```
make serial
./serial -n 256 -i 5000 -s 40
```

**Compile and run**

You can also execute both of the above commands together with default values with `make run_serial`

### Parallel solution
**Compile**

`make parallel_pthreads`
or 
`make parallel_omp`

**Compile and Run**

You can also execute both of the above commands together with default values with `make run_pthreads` or `make run_omp`.

## Visualize
### Plots
`./plot_solution.sh -n [grid_size]`

Plots the program output using [gnuplot](http://gnuplot.sourceforge.net).

Alternatively, you can compile, run, and plot the solution with default values with `make plot_pthreads` or `make plot_omp` .

You can plot the serial solution with `make plot_serial`.

**Example**

`./plot_solution.sh -n 256`

### Video
`make show_pthreads`
or
`make show_omp`

Compiles, runs, and plots the parallel solution with default values and creates a video using [ffmpeg](https://ffmpeg.org).

You can create a video from the serial solution with `make show_serial`.

## Check
`make check_pthreads`
or
`make check_omp`

Compiles and runs the parallel solution with default values and compares the output data to reference data.

You can check the serial solution with `make check_serial`.

## Options
Option | Description | Restrictions | Default value
:------------ | :------------ | :------------ | :------------
**-n** | Number of grid points in one spatial dimension | > 0 | 256
**-i** | Number of iterations | > 0 | 5000
**-s** | Number of iterations between each time the grid state is saved to file | > 0 | 40

## Installing dependencies
**OpenMP and Pthreads**
Your compiler should already include support for OpenMP and pthreads by default.

**gnuplot**

Linux/Ubuntu:

```
sudo apt update
sudo apt install gnuplot
```

MacOSX:

```
brew update
brew install gnuplot
```

**ffmpeg**

Linux/Ubuntu:

```
sudo apt update
sudo apt install ffmpeg
```

MacOSX:

```
brew update
brew install ffmpeg
```
