CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs) -lm

.PHONY: all clean

all: out.txt interp_plot_qspline_deriv.pdf interp_plot_lspline.pdf interp_plot_qspline_int.pdf

interp_plot_qspline_deriv.pdf interp_plot_lspline.pdf interp_plot_qspline_int.pdf: data.txt plot.gpi
	gnuplot plot.gpi

data.txt: main
	./$<

out.txt: main
	./$< > $@


main: main.o linterp.o sinintegral.o qspline.o

main.o linterp.o sinintegral.o: lspline.h

main.o qspline.o: qspline.h

clean: 
	$(RM) main out.txt data.txt *.o *.pdf