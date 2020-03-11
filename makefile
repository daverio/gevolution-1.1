# programming environment
COMPILER     := mpic++
INCLUDE      := -I../LATfield2/
LIB          := -lfftw3 -lm -lhdf5 -lgsl -lgslcblas
# -lcfitsio -lchealpix

# target and source
EXEC         := fRevolution
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5 -DMULTIGRID -DDEBUG_MULTIGRID
# optional compiler settings (LATfield2)
#DLATFIELD2   += -DEXTERNAL_IO # enables I/O server (use with care)
#DLATFIELD2   += -DSINGLE      # switches to single precision, use LIB -lfftw3f

# optional compiler settings (gevolution)
DGEVOLUTION  := -DPHINONLINEAR
DGEVOLUTION  += -DBENCHMARK
DGEVOLUTION  += -DEXACT_OUTPUT_REDSHIFTS
DGEVOLUTION  += -DORIGINALMETRIC
#DGEVOLUTION  += -DVELOCITY      # enables velocity field utilities
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

mapinfo: mapinfo.cpp
	g++ mapinfo.cpp -o mapinfo -std=c++11 -O3

map2fits: map2fits.cpp
	g++ map2fits.cpp -o map2fits -std=c++11 -O3 -I/home/lorerev/Healpix_3.50/include -I/home/lorerev/Healpix_3.50/src/C/subs -I/home/lorerev/cfitsio-3.47/include -L/home/lorerev/Healpix_3.50/src/C/subs -L/home/lorerev/Healpix_3.50/lib -L/home/lorerev/cfitsio-3.47/lib64 -lcfitsio -lchealpix

clean:
	-rm -f $(EXEC) lccat mapinfo map2fits
