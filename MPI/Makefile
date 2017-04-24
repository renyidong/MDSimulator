all: main.c
	mpicc -Wall -O3 -lpthread main.c
debug: main.c
	mpicc -DDEBUG=1 -Wall -g -Og -lpthread  main.c
