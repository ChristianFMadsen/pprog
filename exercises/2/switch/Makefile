.PHONEY:all
all: out.txt
	cat out.txt

out.txt : switch
	./switch > out.txt

switch : switch.c


clean: 
	$(RM) out.txt switch 