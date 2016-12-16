CC = mpic++
CFLAGS  = -g -Wall

default: main
test: init.o parameters.hpp mod_arithm.o test.o
	$(CC) $(CFLAGS) -o test test.o mod_arithm.o init.o
	
main:  init.o  mod_arithm.o sieving.o main.o candidates_main.o candidates_helper.o sparse_binary_vector.o sieving_coord.o sieving_worker.o
	$(CC) $(CFLAGS) -o main main.o init.o mod_arithm.o sieving.o candidates_main.o candidates_helper.o sparse_binary_vector.o sieving_coord.o sieving_worker.o

init.o: init.cpp init.hpp
	$(CC) $(CFLAGS) -c init.cpp

mod_arithm.o: mod_arithm.cpp mod_arithm.hpp 
	$(CC) $(CFLAGS) -c mod_arithm.cpp

sieving_coord.o: sieving_coord.cpp sieving_coord.hpp sieving.hpp
	$(CC) $(CFLAGS) -c sieving_coord.cpp

sieving_worker.o: sieving_worker.cpp 
	$(CC) $(CFLAGS) -c sieving_worker.cpp

sieving.o: sieving.cpp sieving.hpp parameters.hpp  sparse_binary_vector.hpp
	$(CC) $(CFLAGS) -c sieving.cpp	

candidates_main.o: candidates_main.cpp candidates_main.hpp parameters.hpp sparse_binary_vector.hpp sparse_binary_vector.hpp
	$(CC) $(CFLAGS) -c candidates_main.cpp	

candidates_helper.o: candidates_helper.cpp candidates_helper.hpp parameters.hpp sparse_binary_vector.hpp sparse_binary_vector.hpp
	$(CC) $(CFLAGS) -c candidates_helper.cpp	

sparse_binary_vector.o: sparse_binary_vector.cpp sparse_binary_vector.hpp
	$(CC) $(CFLAGS) -c sparse_binary_vector.cpp

test.o: test.cpp mod_arithm.o
	$(CC) $(CFLAGS) -c test.cpp

main.o:  main.cpp init.hpp
	$(CC) $(CFLAGS) -c main.cpp

clean: 
	$(RM) main *.o *~
