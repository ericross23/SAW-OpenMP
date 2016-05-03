CC = gfortran
CFLAGS = -O3 -fopenmp -std='f2003'
TARGETS = project nonreversal dimer

all: $(TARGETS)

project: project.o
	$(CC) $(CFLAGS) -o project project.o

project.o: project.f03
	$(CC) $(CFLAGS) -c project.f03

nonreversal: nonreversal.o
	$(CC) $(CFLAGS) -o nonreversal nonreversal.o

nonreversal.o: nonreversal.f03
	$(CC) $(CFLAGS) -c nonreversal.f03

dimer: dimer.o
	$(CC) $(CFLAGS) -o dimer dimer.o

dimer.o: dimer.f03
	$(CC) $(CFLAGS) -c dimer.f03
