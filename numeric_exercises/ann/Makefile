CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

.PHONY: all

all: annCF.pdf ann2DCF.pdf

annCF.pdf ann2DCF.pdf: out.txt Makefile plot.gpi
	gnuplot plot.gpi

out.txt: main
	./$< > $@ 

main: main.o ann.o minimization.o qr_gs_decomp.o

qr_gs_decomp.o ann.o minimization.o main.o: ann.h

clean:
	rm main *.o *.pdf *.txt