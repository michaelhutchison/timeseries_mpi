# All executables
EXE=player recorder_serial recorder_mpi
# Name of the animation file which will be produced and played back
RECORDFILE=mpi_record.rcd
# Object files required for the playback executable
PLAYEROBJ=graphicslib.o View.o Light.o Scene.o Mouse.o Texture.o vec3.o Model.o 
#  Libraries - Linux
LIBS=-lglut -lGLU -lGL
# Compilers
CPP=g++
MPICC=mpicc
# Size of the world
X=100
Y=100
Z=100
# number of processes to use in MPI implementation
P=4
# Number of frames to produce
F=400
# Number of objects to create in the world
O=2000

# Make and run the MPI recorder
parallel: recorder_mpi 
	rm -f $(RECORDFILE)
	mpirun -np $(P) ./recorder_mpi $(RECORDFILE) $(X) $(Y) $(Z) $(O) $(F)

recorder_mpi: recorder_mpi.o Slice.o Object_mpi.o vec3.o
	$(MPICC) -o recorder_mpi recorder_mpi.o Slice.o Object_mpi.o vec3.o

recorder_mpi.o: recorder_mpi.cpp
	$(MPICC) -c recorder_mpi.cpp

Slice.o: Slice.h Slice.cpp
	$(MPICC) -c Slice.cpp

# Make and run the Serial recorder
serial: recorder_serial 
	rm -f $(RECORDFILE)
	./recorder_serial $(RECORDFILE) $(X) $(Y) $(Z) $(O) $(F)

recorder_serial: recorder_serial.o World.o Object_mpi.o vec3.o
	$(CPP) -o recorder_serial recorder_serial.o World.o Object_mpi.o vec3.o


# make and run the MPI Player
play: player
	./player $(RECORDFILE) $(X) $(Y) $(Z)

#  Generic compile rules
.cpp.o:
	$(CPP) -c -O -Wall $< -IGL

#  Generic compile and link
%: %.cpp graphicslib.a
	$(CPP) -Wall -O3 -o $@ $^ $(LIBS)

#  Delete unwanted files
clean:
	rm -f $(EX) *.o *.a

#  Create archive (include glWindowPos here if you need it)
graphicslib.a: $(PLAYEROBJ) 
	ar -rcs graphicslib.a $^

