CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

.PHONY: all clean 

all: main out.txt plot
	cat out.txt

main:

out.txt: main
	./$< > $@

outAi.txt: main 
	./$<

outBi.txt: main 
	./$<

plot: outAi.txt outBi.txt plot.gpi
	gnuplot plot.gpi

clean: 
	rm main outBi.txt outAi.txt plot.svg out.txt 