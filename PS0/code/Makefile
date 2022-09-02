CC:=gcc
CFLAGS+=
LDLIBS+=-lm

SRC_FILES:=src/shallow_water_serial.c src/argument_utils.c

INC_PATHS:=-Iinc/

.PHONY: clean purge setup run check plot show viz

serial: $(SRC_FILES)
	$(CC) $^ $(CFLAGS) $(LDLIBS) -o $@

clean:
	-rm -f serial 

purge:
	-rm -f serial data/*.bin plots/*.png video/*.mp4

setup:
	-mkdir -p data plots video

run: setup serial
	./serial 

check: setup purge run
	cd check && make && cd ..
	./check/check_solution.sh

plot: setup purge run
	./plot_solution.sh

show: setup purge run viz

viz:
	./plot_solution.sh -n > /dev/null
	ffmpeg -framerate 10 -i plots/%05d.png video/output.mp4 &> /dev/null
	./open_video.sh &> /dev/null




