COMPILER     := CC
INCLUDE      := -I/users/daverio/Documents/LATfield2 -I/users/daverio/local/include
LIB          := -L/users/daverio/local/lib -lfftw3 -lm -lhdf5 -lgsl -lgslcblas -lcfitsio -lchealpix

# target and source
EXEC         := fRevolution
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5 -DMULTIGRID
#-DDEBUG_MULTIGRID
# optional compiler settings (LATfield2)
#DLATFIELD2   += -DEXTERNAL_IO # enables I/O server (use with care)
#DLATFIELD2   += -DSINGLE      # switches to single precision, use LIB -lfftw3f
#DGEVOLUTION  += -DCOLORTERMINAL
#DGEVOLUTION  += -DCHECK_B
#DGEVOLUTION  += -DHAVE_CLASS    # requires LIB -lclass
DGEVOLUTION  += -DHAVE_HEALPIX  # requires LIB -lchealpix

# further compiler options
OPT          := -O3 -std=c++11

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)
	cp $@ /scratch/snx3000/daverio/FR/


lccat: lccat.cpp
	g++ $< -o $@ $(OPT) $(DGEVOLUTION) $(INCLUDE)
