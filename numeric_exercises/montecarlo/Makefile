CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags) 
LDLIBS = $$(gsl-config --libs) -lm
CC = gcc

.PHONY: all
all: out.txt plainMCerror.pdf

out.txt: main Makefile
	./$< > $@

main: montecarlo.o adapInt.o

montecarlo.o adapInt.o: montecarlo.h

errors.txt: main
	./$<

plainMCerror.pdf: plot.gpi errors.txt
	gnuplot plot.gpi

.PHONY:clean
clean:
	$(RM) main *.o *.txt *.pdf *.log
	