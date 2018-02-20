.PHONEY:all
all: out.txt
	cat out.txt

out.txt : hello
	./hello > out.txt

hello : hello.o
	gcc hello.o -o hello

hello.o : hello.c
	gcc -c hello.c -o hello.o

clean: 
	$(RM) out.txt hello hello.o 
