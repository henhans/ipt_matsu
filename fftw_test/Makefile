CC = icpc#g++#mpicxx#icpc
FLAGS = -O3 -g -lfftw3 #-lm #-D_OMP -openmp #-fopenmp
#LFLAGS = #
OBJECTS = mpt.o routines.o main.o


#output.txt: main.exe
#	./main.exe > output.txt

main.exe: $(OBJECTS)
	$(CC) $(FLAGS) $(OBJECTS) -o main.exe

%.o : %.cpp
	$(CC) $(FLAGS) -c -o $@ $< 

clean: 
	rm -f *.o
