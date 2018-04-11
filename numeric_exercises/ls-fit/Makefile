CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags) 
LDLIBS = $$(gsl-config --libs) -lm
CC = gcc


.PHONY: all 
all: leastsquares.pdf

out.txt: main
	./$< > $@


leastsquares.pdf: out.txt Makefile plot.gpi
	gnuplot plot.gpi


main: main.o leastsquares.o qr.o


main.o leastsquares.o qr.o: qr-ls.h


.PHONY: clean
clean:
	$(RM) main *.o *.txt *.pdf 