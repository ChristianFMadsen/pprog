CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags) -O3
LDLIBS += $$(gsl-config --libs)

.PHONY: all
all: annLogistic.pdf annGauss.pdf Makefile

annLogistic.pdf annGauss.pdf: logisticData.txt gaussData.txt plot.gpi
	gnuplot plot.gpi

logisticData.txt gaussData.txt: main
	./$<

main: main.o ann.o minimization.o qr_gs_decomp.o

qr_gs_decomp.o ann.o minimization.o main.o: ann.h

.PHONY: clean
clean:
	rm main *.o *.pdf *.txt