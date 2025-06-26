CC=gcc
ARG=-O3

exec/PInteract: code/PInteract.o
	$(CC) $(ARG) -o exec/PInteract code/PInteract.o -lm
code/PInteract.o: code/PInteract.c
	$(CC) -c $(ARG) code/PInteract.c -o code/PInteract.o
