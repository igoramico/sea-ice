MYFLAGS =-DADAMS_BASHFORTH -DINIT_ICE_FOURIER -DMECHANICAL_CORRELATED 
CFLAGS = -lm -O3 -g
CC = gcc

sea-ice-1d: sea-ice-1d.c sea-ice-1d.h variables.h
	   $(CC) sea-ice-1d.c -o sea-ice-1d $(CFLAGS) $(MYFLAGS)

clean:
	rm -rf sea-ice-1d sea-ice-1d.dSYM

