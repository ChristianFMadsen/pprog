CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONY: all clean

all: main out.txt plot

out.txt: main
	./$< > $@

main:

plot: plot.gpi out.txt
	gnuplot plot.gpi
	
clean: 
	rm main out.txt *.pdf
