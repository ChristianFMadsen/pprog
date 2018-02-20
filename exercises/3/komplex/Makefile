CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
.PHONEY: all clean

all:out.txt 
	cat out.txt

out.txt:main
	./main > out.txt

main:main.o komplex.o

main.o komplex.o : komplex.h

clean:
	$(RM) main main.o komplex.o out.txt
