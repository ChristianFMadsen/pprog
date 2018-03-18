CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONY: all clean

all: main out.txt 
	cat out.txt


main: 


out.txt: main
	./$< > $@

clean: 
	rm main out.txt