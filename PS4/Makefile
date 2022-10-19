CC:=gcc
PARALLEL_CC:=nvcc
CFLAGS+=
LDLIBS+=-lm

SERIAL_SRC_FILES:=src/shallow_water_serial.c src/argument_utils.c

PARALLEL_SRC_FILES:=src/shallow_water_parallel.cu src/argument_utils.cpp

INC_PATHS:=-Iinc/

.PHONY: clean purge setup run check plot show run_serial run_parallel check_serial check_parallel plot_serial plot_parallel show_serial show_parallel viz

parallel: $(PARALLEL_SRC_FILES)
	$(PARALLEL_CC) $^ $(CFLAGS) $(LDLIBS) -o $@

serial: $(SERIAL_SRC_FILES)
	$(CC) $^ $(CFLAGS) $(LDLIBS) -o $@

run: run_parallel
check: check_parallel
plot: plot_parallel
show: show_parallel

clean:
	-rm -f serial parallel 

purge:
	-rm -f serial parallel data/*.bin plots/*.png video/*.mp4

setup:
	-mkdir -p data plots video
	$(MAKE) -C check clean
	$(MAKE) -C check all

run_serial: purge serial
	./serial 

run_parallel: purge parallel
	./parallel

check_serial: purge serial
	./check/check_serial_solution.sh

check_parallel: purge parallel
	./check/check_parallel_solution.sh

plot_serial: purge run_serial
	./plot_solution.sh

plot_parallel: purge run_parallel
	./plot_solution.sh 

show_serial: purge run_serial viz
show_parallel: purge run_parallel viz

viz:
	./plot_solution.sh > /dev/null
	ffmpeg -framerate 10 -i plots/%05d.png video/output.mp4 &> /dev/null
	./open_video.sh &> /dev/null




