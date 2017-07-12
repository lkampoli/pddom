FC = gfortran

FCFLAGS = -Wall -fcheck=all -fopenmp
# FCFLAGS = -fopenmp -O3 -DNDEBUG -march=native -mtune=native -mcmodel=large

LDFLAGS = -L/opt/intel/mkl/lib/intel64

# LIBS = -llapack -lblas -lfftw3 # DESKTOP
LIBS = -lmkl_rt  # LAPTOP
INC = -I/opt/intel/composer_xe_2015.3.187/mkl/include/fftw


# PROGRAMS = benchmark
MODS = special.o quadratures.o discretize.o solver.o halfspace.o


# all: $(PROGRAMS)

benchmark: $(MODS)
benchmark.o: $(MODS)

test: $(MODS)
test.o: $(MODS)



%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.f95
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS) $(LIBS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS) $(LIBS)

%.o: %.f03
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS) $(LIBS)

%.o: %.f08
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS) $(LIBS)


.PHONY: clean cleandata veryclean

clean:
	rm -f *.o *.mod *.MOD
