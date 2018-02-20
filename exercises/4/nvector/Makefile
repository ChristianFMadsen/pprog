CFLAGS = -Wall -std=gnu11
LDLIBS = -lm

.PHONEY: all clean

all: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main: main.o nvector.o

main.o nvector.o: nvector.h

clean: 
	$(RM) main main.o nvector.o out.txt