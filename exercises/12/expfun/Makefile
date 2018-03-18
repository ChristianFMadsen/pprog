CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs) -lm

.PHONY: all clean

all: report.pdf

report.pdf: report.tex plot.tex
	latex $<
	latex $<
	latex $<
	dvipdf report.dvi

plot.tex: out.txt plot.gpi
	gnuplot plot.gpi

out.txt: main
	./$< > $@

main: main.o myexp.o

clean: 
	rm main *.o out.txt plot.tex report.aux report.dvi report.log report.pdf