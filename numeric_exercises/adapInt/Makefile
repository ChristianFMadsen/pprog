CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags) 
LDLIBS = $$(gsl-config --libs) -lm
CC = gcc

.PHONY: all
all: out.txt

out.txt: main Makefile
	./$< > $@

main: adapInt.o

adapInt.o: adapInt.h

.PHONY:clean
clean:
	$(RM) main *.o *.txt