# Select compliler
CMPR = g++

# Create directory
OBJDIR = ./obj/

CPPFLAGS = 

CPP =
CFLAGS = -std=c++17 -Wall  
# For debugging (gdb)
CFLAGS += -g
# For double precision
CFLAGS += -DDOUBLE_PREC
# For OpenMP
# CFLAGS += -fopenmp

OBJ = main.o meshblock.o IC2Dtype.o MUSCL2D.o slopelimiter.o celledges.o riemannS.o savedata.o WENO2D_primvar.o parameters.o diffusivity.o timeIntegral.o
EXEC = main.x
# Multiple static grid stuffs
OBJ += basegrid.o locate_bounds.o boundary.o 
# Making grids adaptive
OBJ += markref.o criteria.o critneighup.o critneighdown.o refine.o coarsen.o updatelpup.o updatelpdown.o admesh.o
# Others
OBJ += misc.o
# Constrained transport for zero solenoidal magnetic field
OBJ += CT2D.o

OBJS = $(addprefix $(OBJDIR), $(OBJ))

$(addprefix ./obj/, %.o): %.cpp
	$(CMPR) -c $(CFLAGS) $(CPPFLAGS) $< -o $@ 

# Link into an executable
main: $(OBJS)
	$(CMPR) $(OBJS) -o $(EXEC)

clean:
	rm -f ./obj/*.o 
	rm -f ./*.x
	rm -f ./data/*.dat 

clean_results:
	rm -f ./data/*.dat ./plots/*.png
