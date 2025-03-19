OBJS = Velocity.o
CC = g++ -O3 -std=c++11
DEBUG = 
CFLAGS = -c $(DEBUG) -I/usr/local/include/blitz
LFLAGS = $(DEBUG) -lfftw3 -lblitz -L/usr/local/lib/



Velocity: $(OBJS)
	$(CC) -o Velocity $(OBJS) $(LFLAGS)

Velocity.o: Velocity.cpp
	$(CC) $(CFLAGS) Velocity.cpp
	

