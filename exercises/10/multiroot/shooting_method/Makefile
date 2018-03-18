CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONY: all clean

all: main out.txt plot.pdf
	cat out.txt


main: main.o Fepsi.o

out.txt data.txt: main
	./$< > $@

plot.pdf: main data.txt 
	gnuplot plot.gpi

clean: 
	rm main *.o *.txt plot.pdf