CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags) 
LDLIBS = $$(gsl-config --libs) -lm
CC = clang


.PHONY: all 
all: out.txt


out.txt: partA partB partC GSLrootfinding Makefile
	./partA > $@
	./partB >> $@
	./GSLrootfinding >> $@
	./partC >> $@


partA: partA.o qr.o rootfinding.o

partB: partB.o qr.o rootfinding.o

partC: partC.o qr.o rootfinding.o

GSLrootfinding: GSLrootfinding.o

partA.o partB.o partC.o qr.o rootfinding.o: qr-rf.h


.PHONY: clean
clean:
	$(RM) main partA partB partC GSLrootfinding *.txt *.o 