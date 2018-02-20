.PHONEY:all
all: out.txt
	cat out.txt

out.txt : exercise-epsilon
	./exercise-epsilon > out.txt

exercise-epsilon : exercise-epsilon.o
	cc exercise-epsilon.o -o exercise-epsilon

exercise-epsilon.o : exercise-epsilon.c
	cc -c exercise-epsilon.c -o exercise-epsilon.o

clean: 
	$(RM) out.txt exercise-epsilon exercise-epsilon.o 
