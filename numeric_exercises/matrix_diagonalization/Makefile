CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags) 
LDLIBS = $$(gsl-config --libs) -lm
CC = clang
matrixDim = 5
numberOfEigVals = 3
matrixDimComparison = 125


.PHONY: all clean

all: outA.txt outB.txt comparisonCyclicVsEigbyeig.txt  matrix_diagonalization_time.pdf


outA.txt: main-cyclic Makefile
	./$< $(matrixDim) > $@

outB.txt: main-eigbyeig Makefile
	./$< $(matrixDim) $(numberOfEigVals) > $@

comparisonCyclicVsEigbyeig.txt: main-cyclic main-eigbyeig-asconly Makefile
	cat /dev/null > comparisonCyclicVsEigbyeig.txt
	@echo "First value is the amount of seconds it takes to diagonalize the entire matrix using the cyclic method, second is how long it takes using the eigbyeig method and the third is how long it takes to find the lowest eigenvalue only:" > comparisonCyclicVsEigbyeig.txt
	\time -f "%U" -ao comparisonCyclicVsEigbyeig.txt ./main-cyclic $(matrixDimComparison) > /dev/null
	\time -f "%U" -ao comparisonCyclicVsEigbyeig.txt ./main-eigbyeig-asconly $(matrixDimComparison) $(matrixDimComparison)  > /dev/null  
	\time -f "%U" -ao comparisonCyclicVsEigbyeig.txt ./main-eigbyeig-asconly $(matrixDimComparison) 1  > /dev/null 


time.txt: main-cyclic Makefile 
	cat /dev/null > time.txt #Cleans time.txt
	for number in `seq 50 5 225`; do \time --format "$$number \t %U" --append --output $@ ./$< $$number > /dev/null; done

main-cyclic: jacobi.o main-cyclic.o printMatrix.o

main-eigbyeig: jacobi.o main-eigbyeig.o printMatrix.o

main-eigbyeig-asconly: jacobi.o main-eigbyeig-asconly.o printMatrix.o 

main-cyclic.o main-eigbyeig.o jacobi.o printMatrix.o: md.h

matrix_diagonalization_time.pdf: plot.gpi time.txt
	gnuplot plot.gpi

clean:
	$(RM) main-cyclic main-eigbyeig main-eigbyeig-asconly *.o *.txt *.pdf *.log





