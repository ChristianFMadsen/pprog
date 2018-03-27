CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags) 
LDLIBS = $$(gsl-config --libs) -lm
CC = clang

.PHONY: all clean

all: out.txt

out.txt: main
	./$< > $@

main: main.o qr_gs_decomp.o

main.o qr_gs_decomp.o: lin_eq.h

clean: 
	$(RM) main out.txt *.o