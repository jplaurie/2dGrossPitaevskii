OBJS = sound.o
#this is for linux machines
CC = g++ -Wall -O3 -std=c++14 -fopenmp
CFLAGS = -c -I/usr/local/include
LFLAGS = $(DEBUG) -lfftw3 -lfftw3_threads -lopenblas -L/usr/local/lib -lm -lpthread -llapack -larmadillo


#this is for aston math cluster
#CC = g++ -Wall -O3 -std=c++14 -fopenmp
#CFLAGS = -c -I/apps/lib64/openblas/include -I/home/lauriej/libs/include -I/home/lauriej/libs/    usr/local/include
#LFLAGS = -lfftw3 -lfftw3_threads -L/appa/lib64/openbals/lib -L/home/lauriej/libs/lib -L/home/    lauriej/libs/usr/local/lib64 -lm -lpthread -lopenblas -llapack -larmadillo



#The following is for mac os x
#CC = g++-mp-10 -Wall -O3 -std=c++14 -fopenmp
#DEBUG = 
#CFLAGS = -c $(DEBUG) -I/opt/local/include
#LFLAGS = $(DEBUG) -L/opt/local/lib -framework Accelerate -march=native -lfftw3 -lfftw3_threads -lopenblas -lm -lpthread -larmadillo



sound: $(OBJS)
	$(CC) -o sound $(OBJS) $(LFLAGS)

sound.o: sound.cpp
	$(CC) $(CFLAGS) sound.cpp
	

