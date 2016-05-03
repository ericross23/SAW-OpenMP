CC = gfortran
CFLAGS = -O3 -fopenmp -std='f2003'

project: project.o
	$(CC) $(CFLAGS) -o project project.o

project.o: project.f03
	$(CC) $(CFLAGS) -c project.f03

serial: serial.o
	$(CC) $(CFLAGS) -o serial serial.o

serial.o: serial.f03
	$(CC) $(CFLAGS) -c serial.f03
