OBJS = main.o linearoperator.o nonlinearoperator.o read.o energy.o intensity.o spectrum.o forcing.o flux.o
#this is for linux machines
CC = g++ -Wall -O3 -std=c++14 -fopenmp

CUFFTW_INC_PATH = -I/usr/local/cuda-11.7/targets/x86_64-linux/include
CUFFTW_LIB_PATH = -L/usr/local/cuda-11.7/targets/x86_64-linux/lib

DEBUG = 
CFLAGS = -c $(DEBUG) 
LFLAGS = $(DEBUG) -lcufftw -lcufft -lopenblas -lm -lpthread -llapack -larmadillo $(CUFFTW_INC_PATH) $(CUFFTW_LIB_PATH)







quantum2d: $(OBJS)
	$(CC) -o quantum2d $(OBJS) $(LFLAGS)

main.o: main.cpp const.h functions.h
	$(CC) $(CFLAGS) $(LFLAGS) main.cpp
	
nonlinearoperator.o: nonlinearoperator.cpp const.h functions.h
	$(CC) $(CFLAGS) $(LFLAGS) nonlinearoperator.cpp

linearoperator.o: linearoperator.cpp const.h functions.h
	$(CC) $(CFLAGS) $(LFLAGS) linearoperator.cpp

read.o: read.cpp const.h functions.h
	$(CC) $(CFLAGS) $(LFLAGS) read.cpp
    
energy.o: energy.cpp const.h functions.h
	$(CC) $(CFLAGS) $(LFLAGS) energy.cpp

intensity.o: intensity.cpp const.h functions.h
	$(CC) $(CFLAGS) $(LFLAGS) intensity.cpp

spectrum.o: spectrum.cpp const.h functions.h
	$(CC) $(CFLAGS) $(LFLAGS) spectrum.cpp

forcing.o: forcing.cpp const.h functions.h
	$(CC) $(CFLAGS) $(LFLAGS) forcing.cpp

#Filter.o: Filter.cpp Const.h Functions.h
#	$(CC) $(CFLAGS) Filter.cpp

flux.o: flux.cpp const.h functions.h
	$(CC) $(CFLAGS) $(LFLAGS) flux.cpp

clean:
	rm -rfv *.o ./Output/ *~ quantum2d

