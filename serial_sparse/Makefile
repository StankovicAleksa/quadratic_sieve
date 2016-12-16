CC = g++
CFLAGS  = -g -Wall

default: main
test: init.o parameters.hpp mod_arithm.o test.o
	$(CC) $(CFLAGS) -o test test.o mod_arithm.o init.o
	
main:  init.o  mod_arithm.o sieving.o main.o candidates.o sparse_binary_vector.o
	$(CC) $(CFLAGS) -o main main.o init.o mod_arithm.o sieving.o candidates.o sparse_binary_vector.o

init.o: init.cpp init.hpp
	$(CC) $(CFLAGS) -c init.cpp

mod_arithm.o: mod_arithm.cpp mod_arithm.hpp 
	$(CC) $(CFLAGS) -c mod_arithm.cpp

sieving.o: sieving.cpp sieving.hpp parameters.hpp  sparse_binary_vector.hpp
	$(CC) $(CFLAGS) -c sieving.cpp	

candidates.o: candidates.cpp candidates.hpp parameters.hpp sparse_binary_vector.hpp sparse_binary_vector.hpp
	$(CC) $(CFLAGS) -c candidates.cpp	

sparse_binary_vector.o: sparse_binary_vector.cpp sparse_binary_vector.hpp
	$(CC) $(CFLAGS) -c sparse_binary_vector.cpp

test.o: test.cpp mod_arithm.o
	$(CC) $(CFLAGS) -c test.cpp

main.o:  main.cpp init.hpp
	$(CC) $(CFLAGS) -c main.cpp

clean: 
	$(RM) main *.o *~
