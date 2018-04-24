CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags) 
LDLIBS = $$(gsl-config --libs) -lm
CC = gcc

.PHONY: all
all: out.txt sinPlot.pdf

out.txt: main Makefile
	./$< > $@

main: ode.o

ode.o: ode.h

sin.txt: main
	./$<

sinPlot.pdf: plot.gpi sin.txt
	gnuplot plot.gpi

.PHONY:clean
clean:
	$(RM) main *.o *.txt *.pdf