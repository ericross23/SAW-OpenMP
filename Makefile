CC = gfortran
CFLAGS = -O3

project: project.o
	$(CC) -o project project.o

project.o: project.f95
	$(CC) -c project.f95
