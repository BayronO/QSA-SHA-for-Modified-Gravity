FC      = gfortran
MOD     = mod
INCDIR = .
#INCDIR = /usr/include/gsl
INCDIR1 = .
#INCDIR1 = /cygdrive/c/Users/ragnar/Desktop/fgsl-1.1.0
LIBDIR = .
#LIBDIR = /cygdrive/c/Users/ragnar/Desktop/fgsl-1.1.0/.libs
OBJS    = dverk.o 
FFLAGS  = -O3 -ffree-line-length-5000 -I$(INCDIR) -I$(INCDIR1) -fno-second-underscore -fPIC -g
#LDFLAGS = -L$(LIBDIR) -lfgsl -lgsl 
LDFLAGS = -L$(LIBDIR)

.SUFFIXES: .f90 

.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@

default: designer

designer: $(OBJS) designer.o
	$(FC) $(FFLAGS) -o $@ $@.o $(OBJS) $(LDFLAGS)

designer_bat: $(OBJS) designer_bat.o
	$(FC) $(FFLAGS) -o $@ $@.o $(OBJS) $(LDFLAGS)
	
HS: $(OBJS) HS.o
	$(FC) $(FFLAGS) -o $@ $@.o $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o *~

tidy: clean
	rm -f sample
	rm -f *.$(MOD)
