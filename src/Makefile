OBJS = main.o linear.o nonlinear.o read.o compute.o forcing.o print.o
#this is for linux machines
#CC = g++ -Wall -O3 -std=c++14 -fopenmp
#CFLAGS = -c -I/usr/local/include
#LFLAGS = -lfftw3 -lfftw3_threads -lopenblas -L/usr/local/lib -lm -lpthread -llapack -larmadillo


#The following is for mac os x
CC = g++-mp-14 -O3 -std=c++23 -fopenmp
DEBUG = -Wall 
CFLAGS = -c $(DEBUG) -I/opt/local/include
LFLAGS = $(DEBUG) -L/opt/local/lib -framework Accelerate -march=native -lfftw3 -lfftw3_threads -lopenblas -lm -lpthread -larmadillo



gp2d: $(OBJS)
	$(CC) -o gp2d $(OBJS) $(LFLAGS)

main.o: main.cpp const.h linear.h read.h print.h compute.h forcing.h
	$(CC) $(CFLAGS) main.cpp
	
nonlinear.o: nonlinear.cpp const.h 
	$(CC) $(CFLAGS) nonlinear.cpp

linear.o: linear.cpp const.h 
	$(CC) $(CFLAGS) linear.cpp

read.o: read.cpp const.h 
	$(CC) $(CFLAGS) read.cpp
    
compute.o: compute.cpp const.h
	$(CC) $(CFLAGS) compute.cpp

print.o: print.cpp const.h
	$(CC) $(CFLAGS) print.cpp

forcing.o: forcing.cpp const.h 
	$(CC) $(CFLAGS) forcing.cpp

#Filter.o: Filter.cpp Const.h Functions.h
#	$(CC) $(CFLAGS) Filter.cpp

clean:
	rm -rfv *.o ./parameters.txt ./output/ *~ gp2d

