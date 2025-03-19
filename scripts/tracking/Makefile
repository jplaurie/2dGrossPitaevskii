OBJS = vortex_tracking.o


#for linux
CC = g++ -Wall -O3 -std=c++14 -fopenmp
CFLAGS = -c -I/usr/local/include
LFLAGS = -lfftw3 -lfftw3_threads -lopenblas -L/usr/local/lib -lm -lpthread -llapack -larmadillo



#for mac
#CC = g++-mp-8 -O2 -std=c++14 -fopenmp
#DEBUG = 
#CFLAGS = -c $(DEBUG) -I/opt/local/include
#LFLAGS = $(DEBUG) -L/opt/local/lib -framework Accelerate -march=native -lfftw3 -lfftw3_threads -lopenblas -lm -lpthread -larmadillo



vortex_tracking: $(OBJS)
	$(CC) -o vortex_tracking $(OBJS) $(LFLAGS)

vortex_tracking.o: vortex_tracking.cpp
	$(CC) $(CFLAGS) vortex_tracking.cpp
	

