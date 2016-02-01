%.o: %.c
	$(CC) $(CCOMPILEFLAGS) $*.c
# compiler
CC = icc
# CC = icc
CCOMPILEFLAGS  =  -c -O3 
LIBS    = 
INCS    = 
LINK_SUFFIX  =     -lm

### objects
OBJ = generic.o gridmaps.o gridmanip.o io_fieldview.o  io_ugrid.o  VTK_Out.o flowsolver.o main.o
###

### compile and link
EXE = ns3d

all: $(EXE)
ns3d:  $(OBJ)
	$(CC)  -o $@ $(OBJ) $(LINK_SUFFIX)
clean:
	rm ns3d *.o *.vtk *~ *.dat
