CFLAGS = -Wall -Werror -std=gnu11
LDLIBS = -lm

obj1 = main-stdio.o
obj2 = main-cmdline.o


.PHONEY: clean

main-stdio: $(obj1)

main-cmdline: $(obj2) 

out.io.txt: main-stdio input.data.txt
	cat input.data.txt | ./main-stdio > out.io.txt

out.cmd.txt: main-cmdline input.data.txt
	./main-cmdline `cat input.data.txt` > out.cmd.txt

plot: plot.gpi
	gnuplot plot.gpi

clean: 
	$(RM) $(obj1) $(obj2) main-stdio main-cmdline out.io.txt out.cmd.txt plot.svg
