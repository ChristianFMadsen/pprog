CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags) 
LDLIBS = $$(gsl-config --libs) -lm
CC = clang


.PHONY: all 
all: out.txt


out.txt: partA partB GSLrootfinding Makefile
	./partA > $@
	./partB >> $@
	./GSLrootfinding >> $@


partA: partA.o qr.o rootfinding.o


partB: partB.o qr.o rootfinding.o

GSLrootfinding: GSLrootfinding.o

partA.o partB.o qr.o rootfinding.o: qr-rf.h


.PHONY: clean
clean:
	$(RM) main partA partB GSLrootfinding *.txt *.o 