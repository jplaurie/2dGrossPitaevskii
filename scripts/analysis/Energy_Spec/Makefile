OBJS = Energy.o 
CC = g++ -O2 -std=c++11
DEBUG = 
CFLAGS = -c $(DEBUG) -I/usr/local/include/blitz
LFLAGS = $(DEBUG) -lfftw3 -lblitz -L/usr/local/lib/



Energy: $(OBJS)
	$(CC) -o Energy $(OBJS) $(LFLAGS)

Energy.o: Energy.cpp 
	$(CC) $(CFLAGS) Energy.cpp
	
clean:
	rm -rfv *.o*~ Energy ./kin*

