# Name of the animation file which will be produced and played back
DATA=mpi_record.rcd
#  Libraries - Linux
LIBS=-lglut -lGLU -lGL
# Compilers
CPP=g++
MPICC=mpicc

# make and run the MPI recorder
record: recorder_mpi.cpp Slice.h Slice.cpp Object_mpi.h Object_mpi.cpp vec3.h vec3.cpp 
	rm -f $(DATA)
	$(MPICC) -c recorder_mpi.cpp
	$(MPICC) -c Slice.cpp
	$(MPICC) -c Object_mpi.cpp
	$(MPICC) -c vec3.cpp
	$(MPICC) -o recorder_mpi recorder_mpi.o Slice.o Object_mpi.o vec3.o
	mpirun -np 4 ./recorder_mpi $(DATA) 

# make and run the MPI Player
play: player

#  Generic compile rules
.cpp.o:
	$(CPP) -c -O -Wall $< -IGL

#  Generic compile and link
%: %.cpp graphicslib.a
	$(CPP) -Wall -O3 -o $@ $^ $(LIBS)

#  Delete unwanted files
clean:
	rm -f $(EX) $(DATA) *.o *.a

#  Create archive (include glWindowPos here if you need it)
graphicslib.a: $(OBJ) 
	ar -rcs graphicslib.a $^

