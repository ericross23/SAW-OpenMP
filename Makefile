CC = gfortran
CFLAGS = -O3 -fopenmp

project: project.o
	$(CC) $(CFLAGS) -o project project.o

project.o: project.f95
	$(CC) $(CFLAGS) -c project.f95
