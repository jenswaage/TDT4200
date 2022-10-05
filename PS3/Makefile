CC:=gcc
PARALLEL_CC:=gcc
CFLAGS+= -fopenmp
LDLIBS+=-lm -pthread

SERIAL_SRC_FILES:=src/shallow_water_serial.c src/argument_utils.c

OMP_SRC_FILES:=src/shallow_water_parallel_omp.c src/argument_utils.c
PTHREADS_SRC_FILES:=src/shallow_water_parallel_pthreads.c src/argument_utils.c

INC_PATHS:=-Iinc/

.PHONY: clean purge setup run check plot show run_serial run_parallel check_serial check_parallel plot_serial plot_parallel show_serial show_parallel viz

pthreads: $(PTHREADS_SRC_FILES)
	$(PARALLEL_CC) $^ $(CFLAGS) $(LDLIBS) -o $@

omp: $(OMP_SRC_FILES)
	$(PARALLEL_CC) $^ $(CFLAGS) $(LDLIBS) -o $@

serial: $(SERIAL_SRC_FILES)
	$(CC) $^ $(CFLAGS) $(LDLIBS) -o $@

run: run_pthreads
check: check_pthreads
plot: plot_pthreads
show: show_pthreads

clean:
	-rm -f serial pthreads omp 

purge:
	-rm -f serial pthreads omp data/*.bin plots/*.png video/*.mp4

setup:
	-mkdir -p data plots video
	$(MAKE) -C check clean
	$(MAKE) -C check all

run_serial: purge serial
	./serial

run_omp: purge omp
	./omp

run_pthreads: purge pthreads
	./pthreads

check_serial: purge serial
	./check/check_serial_solution.sh

check_pthreads: purge pthreads
	./check/check_pthreads_solution.sh

check_omp: purge omp
	./check/check_omp_solution.sh

plot_serial: purge run_serial
	./plot_solution.sh

plot_pthreads: purge run_pthreads
	./plot_solution.sh

plot_omp: purge run_omp
	./plot_solution.sh

show_serial: purge run_serial viz
show_pthreads: purge run_pthreads viz
show_omp: purge run_omp viz

viz:
	./plot_solution.sh > /dev/null
	ffmpeg -framerate 10 -i plots/%05d.png video/output.mp4 &> /dev/null
	#./open_video.sh &> /dev/null




