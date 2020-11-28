# Makefile for raysum programs. Optimize as much as possible
FF = gfortran
FFLAGS = -O3 -w
LDFLAGS = -O3 -w

# Subroutines
EIG = eigenvec.o matrixops.o eispack-cg.o
RAYSUM = raysum.o $(EIG)
IO = readwrite.o
TRACE = phaselist.o buildmodel.o trace.o

# Implicit rule for Fortran compilation. Changing the parameter header
# file recompiles everything.
.f.o: params.h
	$(FF) $(FFLAGS) -c $<

default: all

all: seis-spread

seis-spread: seis-spread.o $(RAYSUM) $(IO) $(TRACE)
	$(FF) $(LDFLAGS) -o seis-spread seis-spread.o $(RAYSUM) $(IO) $(TRACE)

clean:
	/bin/rm -f *.o *.core seis-spread

tidy: all clean