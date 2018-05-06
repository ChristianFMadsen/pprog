CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)


annCF.pdf: out.txt Makefile plot.gpi
	gnuplot plot.gpi

out.txt: main
	./$< > $@ 

main: main.o ann.o minimization.o qr_gs_decomp.o

qr_gs_decomp.o ann.o minimization.o main.o: ann.h

clean:
	rm main *.o *.pdf *.txt