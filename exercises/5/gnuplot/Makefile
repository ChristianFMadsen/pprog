CFLAGS = -Wall -Werror -std=gnu11
LDLIBS = -lm



.PHONEY: clean all

all: main-stdio main-cmdline out.io.txt out.cmd.txt input.data.txt plot 

main-stdio:

main-cmdline:

input.data.txt: 
	seq 0 0.02 5 > $@

out.io.txt: main-stdio input.data.txt
	cat input.data.txt | ./main-stdio > out.io.txt

out.cmd.txt: main-cmdline input.data.txt
	./main-cmdline `cat input.data.txt` > out.cmd.txt

plot: plot.gpi input.data.txt
	gnuplot plot.gpi

clean: 
	$(RM)  main-stdio main-cmdline out.io.txt out.cmd.txt plot.svg main-stdio.o main-cmdline.o input.data.txt
