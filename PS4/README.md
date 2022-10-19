# TDT4200 Problem set 4: Intro to CUDA 

## Finite difference approximation of the 1D shallow water equations using CUDA 
In this assignment you will take a serial implementation of the Finite Difference Method (FDM) for solving the 1D shallow water equations and write a GPU version using CUDA. Information on parallelization with CUDA is described in the lecture slides.

The serial solution can be found in `shallow_water_serial.c` and should be kept as a reference. A copy of the serial solution can be found in `shallow_water_parallel.cu`, in which you should write your parallel implementation.

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
./serial -n 1024 -i 100000 -s 1000
```

**Compile and run**

You can also execute both of the above commands together with default values with `make run_serial`

### Parallel solution
**Compile**

`make parallel`

**Compile and Run**

You can also execute both of the above commands together with default values with `make run`.

## Visualize
### Plots
`./plot_solution.sh -n [grid_size]`

Plots the program output using [gnuplot](http://gnuplot.sourceforge.net).

Alternatively, you can compile, run, and plot the solution with default values with `make plot` .

You can plot the serial solution with `make plot_serial`.

**Example**

`./plot_solution.sh -n 1024`

### Video
`make show`

Compiles, runs, and plots the parallel solution with default values and creates a video using [ffmpeg](https://ffmpeg.org).

You can create a video from the serial solution with `make show_serial`.

## Check
`make check`

Compiles and runs the parallel solution with default values and compares the output data to reference data.

You can check the serial solution with `make check_serial`.

## Options
Option | Description | Restrictions | Default value
:------------ | :------------ | :------------ | :------------ 
**-n** | Number of grid points in one spatial dimension | > 0 | 1024
**-i** | Number of iterations | > 0 | 100000
**-s** | Number of iterations between each time the grid state is saved to file | > 0 | 1000
**-np**| Number of processes (MPI option) | > 0 | 4

## Installing dependencies
If you do not have an Nvidia GPU available to use, or do not have the correct environment, you can use our servers through Snotra. See the Sources & Syllabus post on how to use the cluster.


**CUDA**

Linux/Ubuntu:

```
sudo apt update
sudo apt install -y nvidia-cuda-toolkit
```

**gnuplot**

Linux/Ubuntu:

```
sudo apt update
sudo apt install gnuplot
```

**ffmpeg**

Linux/Ubuntu:

```
sudo apt update
sudo apt install ffmpeg
```
