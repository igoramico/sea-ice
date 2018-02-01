MYFLAGS =-DDISJOINING 
CFLAGS = -lm -O3 -g
CC = gcc

sea-ice: sea-ice.c sea-ice.h variables.h
	   $(CC) sea-ice.c -o sea-ice $(CFLAGS) $(MYFLAGS)

clean:
	rm -f sea-ice sea-ice.dSYM

