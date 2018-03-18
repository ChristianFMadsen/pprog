CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONY: all clean

all: report.pdf

report.pdf: plot.tex report.tex
	latex report.tex
	latex report.tex
	latex report.tex
	dvipdf report.dvi

plot.tex: plot.gpi data.txt
	gnuplot plot.gpi

data.txt: main Makefile
	./$< -3 3 0.01 > $@

 

clean: 
	rm main data.txt plot.tex report.dvi report.log report.pdf report.aux
