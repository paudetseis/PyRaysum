# Makefile for raysum programs.
FF = gfortran
FFLAGS = -C -g
LDFLAGS = -C -g

EIG = eigenvec.o matrixops.o eispack-cg.o
RAYSUM = raysum.o $(EIG)
IO = readwrite.o
TRACE = phaselist.o buildmodel.o trace.o
MISFIT = misfit.o

# implicit rule for Fortran compilation. Changing the parameter header
# file recompiles everything.
.f.o: params.h
	$(FF) $(FFLAGS) -c $<

default: all

all: seis-spread

seis-spread: seis-spread.o $(RAYSUM) $(IO) $(TRACE)
	$(FF) $(LDFLAGS) -o seis-spread seis-spread.o $(RAYSUM) \
                             $(IO) $(TRACE)
	cp seis-spread ~/opt/anaconda3/envs/prs/bin

clean:
	/bin/rm -f *.o *.core seis-spread

