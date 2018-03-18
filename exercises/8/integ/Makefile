CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONY: all clean

all: main out.txt plot
	cat out.txt

main:


out.txt: main
	./$< > $@

energies.txt: main
	./$< > $@	

plot: plot.gpi energies.txt
	gnuplot plot.gpi

clean: 
	rm main energies.txt out.txt plot