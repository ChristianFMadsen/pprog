.PHONEY:all
all: out.txt
	cat out.txt

out.txt : math
	./math > out.txt

math : math.o 
	gcc math.o -o math

math.o : math.c
	gcc -c math.c -o math.o

clean: 
	$(RM) out.txt math math.o 