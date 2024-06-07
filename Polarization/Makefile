OBJS = polar.o 
CC = g++ -O2 -std=c++11
DEBUG = 
CFLAGS = -c $(DEBUG) -I/usr/local/include/blitz
LFLAGS = $(DEBUG) -lfftw3 -lblitz -L/usr/local/lib/



polar: $(OBJS)
	$(CC) -o polar $(OBJS) $(LFLAGS)

polar.o: polar.cpp 
	$(CC) $(CFLAGS) polar.cpp
	
clean:
	rm -rfv *.o*~ polar

