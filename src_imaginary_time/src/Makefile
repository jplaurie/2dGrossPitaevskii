OBJS = main.o linearoperator.o nonlinearoperator.o read.o energy.o intensity.o spectrum.o
#this is for linux machines
CC = g++ -Wall -O3 -std=c++14 -fopenmp
DEBUG = 
CFLAGS = -c $(DEBUG) -I/usr/local/include/
LFLAGS = $(DEBUG) -lfftw3 -lfftw3_threads -larmadillo -L/usr/local/lib/ -lm -lpthread

#this is for aston math cluster
#CC = g++ -Wall -O3 -std=c++14 -fopenmp
#CFLAGS = -c -I/apps/lib64/openblas/include -I/home/lauriej/libs/include -I/home/lauriej/libs/    usr/local/include
#LFLAGS = -lfftw3 -lfftw3_threads -L/appa/lib64/openbals/lib -L/home/lauriej/libs/lib -L/home/    lauriej/libs/usr/local/lib64 -lm -lpthread -lopenblas -llapack -larmadillo



#The following is for mac os x
#CC = g++-mp-8 -Wall -O2 -std=c++14 -fopenmp
#DEBUG = 
#CFLAGS = -c $(DEBUG) -I/opt/local/include
#LFLAGS = $(DEBUG) -L/opt/local/lib -framework Accelerate -march=native -lfftw3 -lfftw3_threads -lopenblas -lm -lpthread -larmadillo



imaginary_time_gp: $(OBJS)
	$(CC) -o imaginary_time_gp $(OBJS) $(LFLAGS)

main.o: main.cpp const.h functions.h
	$(CC) $(CFLAGS) main.cpp
	
nonlinearoperator.o: nonlinearoperator.cpp const.h functions.h
	$(CC) $(CFLAGS) nonlinearoperator.cpp

linearoperator.o: linearoperator.cpp const.h functions.h
	$(CC) $(CFLAGS) linearoperator.cpp

read.o: read.cpp const.h functions.h
	$(CC) $(CFLAGS) read.cpp
    
energy.o: energy.cpp const.h functions.h
	$(CC) $(CFLAGS) energy.cpp

intensity.o: intensity.cpp const.h functions.h
	$(CC) $(CFLAGS) intensity.cpp

spectrum.o: spectrum.cpp const.h functions.h
	$(CC) $(CFLAGS) spectrum.cpp

clean:
	rm -rfv *.o ./Output/ *~ quantum2d

