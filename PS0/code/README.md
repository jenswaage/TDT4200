# TDT4200 Problem set 0: Intro to C

## Finite difference approximation of the 1D shallow water equations

In this assignment you will write a serial implementation of the Finite Difference Method (FDM) for solving the 1D shallow water equations. This serial implementation will be the baseline for future problem sets. Information on solving the shallow water equations using FDM is described in the lecture slides.

The skeleton for your serial implementation can be found in `shallow_water_serial.c`.

## Run

### Setup

`make setup`

Creates folders `data`, `plots` and `video`.

- `data`: contains output from the simulation
- `plots`: contains output from plotting
- `video`: contains output from video generation

### Compile

`make serial`

### Run
`./serial -n [grid_size] -i [max_iteration] -s [snapshot_frequency]`

### Example
```
make serial
./serial -n 1024 -i 100000 -s 1000
```

### Compile and run
You can also execute both of the above commands together with default values with `make run`.

## Visualize
### Plots
`./plot_solution.sh -n [grid_size]`

Plots the program output using [gnuplot](http://gnuplot.sourceforge.net).

Alternatively, you can compile, run, and plot the solution with default values with `make plot` .

**Example**

`./plot_solution.sh -n 1024`

### Video
`make show`

Compiles, runs, and plots the solution with default values and creates a video using [ffmpeg](https://ffmpeg.org).

## Check
`make check`

Compiles and runs the solution with default values and compares the output data to reference data.

## Options
Option | Description | Restrictions | Default value
:------------ | :------------ | :------------ | :------------ 
**-n** | Number of grid points in one spatial dimension | > 0 | 1024
**-i** | Number of iterations | > 0 | 100000
**-s** | Number of iterations between each time the grid state is saved to file | > 0 | 1000

## Installing dependencies
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
