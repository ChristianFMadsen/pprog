CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags) 
LDLIBS = $$(gsl-config --libs) -lm
CC = clang


.PHONY: all 
all: out.txt plot.pdf

out.txt: partA partB Makefile
	./$< > $@
	./partB >> $@

data.txt: partB
	./$<

plot.pdf: plot.gpi data.txt Makefile 
	gnuplot plot.gpi

partA: minimization.o qr_gs_decomp.o partA.o 

partB: minimization.o qr_gs_decomp.o partB.o 

partA.o partB.o minimization.o qr_gs_decomp.o main.o: minimization.h

.PHONY: clean
clean:
	$(RM) *.o *.pdf *.txt partA partB