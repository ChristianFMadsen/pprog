.PHONEY:all
all: out.txt
	cat out.txt

out.txt : exercise-precision
	./exercise-precision > out.txt

exercise-precision : exercise-precision.o function-precision.o
	cc exercise-precision.o function-precision.o -o exercise-precision

exercise-precision.o : exercise-precision.c
	cc -c exercise-precision.c -o exercise-precision.o

function-precision.o : function-precision.c
	cc -c function-precision.c -o function-precision.o

clean: 
	$(RM) out.txt exercise-precision.o function-precision.o exercise-precision