all: mex libspartan.a 618.a demo sfddemo

# These BASE variables should be altered to point to the AMD, UFCONFIG 
# and UMFPACK paths
AMD_BASE       = /home/chris/solvers/SuiteSparse/AMD
UFCONFIG_BASE  = /home/chris/solvers/SuiteSparse/UFconfig/
UMFPACK_BASE   = /home/chris/solvers/SuiteSparse/UMFPACK

CC             = gcc
CFLAGS         = -O3 -fPIC
MEX_EXTENSION  = mexa64
FC             = g95
FFLAGS         = -O3 -fPIC

# Should not need to edit past this point
AMD_INC        = -I$(AMD_BASE)/Include
UFCONFIG_INC   = -I$(UFCONFIG_BASE)
UMFPACK_INC    = -I$(UMFPACK_BASE)/Include $(AMD_INC) $(UFCONFIG_INC)
UMFPACK_LIBDIR = -L$(UMFPACK_BASE)/Lib 
UMFPACK_LIBS   = -lumfpack -lamd -lblas

SPARTAN_LIB    = -L. -Wl,-Bstatic -lspartan -Wl,-Bdynamic
MEXDOGLEGFILES =  spartan_mex.c

LIBSPARTANFILES = spartan_blas.c \
	          spartan_malloc.c \
		  spartan_print.c \
	          spartan_dogleg.c \
	          spartan_linsolve.c \
	          spartan_gatxpy.c \
	          spartan_gaxpy.c \
                  spartan_norms.c \
	          spartan_util.c \
	          spartan_entry.c \
                  spartan_compress.c \
	          spartan_cumsum.c \
	          spartan_solve.c \
		  spartan_sfdjac.c

MEX_FILES = spartan_dogleg.$(MEX_EXTENSION)

mex: $(MEX_FILES)

$(MEX_FILES) : $(MEXDOGLEGFILES) $(LIBSPARTANFILES)
	#Note the LP64 symbol needs to be definied if -largeArrayDims is used
	mex -v -largeArrayDims -g3 -DLP64 -output spartan_dogleg $(MEXDOGLEGFILES) \
	$(LIBSPARTANFILES) $(UMFPACK_INC) $(UMFPACK_LIBDIR) $(UMFPACK_LIBS) 


libspartan.a: $(LIBSPARTANFILES)
	$(CC) $(CFLAGS) -c $(LIBSPARTANFILES) $(UMFPACK_INC) 
	ar -rs $@ $(LIBSPARTANFILES:.c=.o)

demo: demo.c libspartan.a
	$(CC) $(CFLAGS) demo.c -o $@ $(UMFPACK_INC) \
	$(SPARTAN_LIB) $(UMFPACK_LIBDIR) $(UMFPACK_LIBS) 

sfddemo: fddemo.c 618.f spartan_columncolor.c libspartan.a
	$(CC) $(CFLAGS) -c spartan_columncolor.c fddemo.c
	$(FC) $(FFLAGS) -c 618.f 
	$(CC) $(CFLAGS) 618.o spartan_columncolor.o fddemo.o -o $@ \
	$(SPARTAN_LIB) $(UMFPACK_LIBDIR) $(UMFPACK_LIBS) 

618.a: 618.f 
	$(FC) $(FFLAGS) -c 618.f
	ar -rs $@ 618.o


run: demo
	./demo

clean:
	-rm *.o	
	-rm demo
	-rm sfddemo
	-rm $(MEX_FILES)
	-rm libspartan.a