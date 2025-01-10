PS_PATH = ../FDPS/src/
INC = -I$(PS_PATH)
CC = time g++
CFLAGS = -std=c++20 -O3 -Wall

#test_kepler.out: test_kepler.cpp
#	$(CC) $(CFLAGS) $(INC) $< -o $@

f_Euler.out: f_Euler.cpp
	$(CC) $(CFLAGS) $(INC) $< -o $@
