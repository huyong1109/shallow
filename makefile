# -------------------Compiler & Flags------------------------------------------
F90COMPILER = ifort
F90LINKER   = ifort 
F90COMPILEOPTS=-O2 -fPIC -m64 -convert big_endian
F90COMPILEPATHS=-I/home/share/cesm/software/esm-soft/include
F90LINKOPTS= -m64 -Wl,--no-as-needed 
F90LINKPATHS=-L/home/share/cesm/software/esm-soft/lib 
F90LINKLIBS=-lnetcdff -lnetcdf -lpnetcdf -lhdf5 -lhdf5_hl -ldl -lm -lz -lcurl -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

# -----------------------------------------------------------------------------
SRCDIR =  src

VPATH= $(SRCDIR)

.SUFFIXES: .F .f .f90 .F90 

%.o : %.f90
	$(F90COMPILER) -c $(F90COMPILEOPTS) $(F90COMPILEPATHS) $<

%.o : %.F90
	$(F90COMPILER) -c $(F90COMPILEOPTS) $(F90COMPILEPATHS) $(F90COMPILECPPFLAGS) $<
 
%.o : %.f
	$(F90COMPILER) -c $(F90COMPILEOPTS) $(F90COMPILEPATHS) $<
 
%.o : %.F
	$(F90COMPILER) -c $(F90COMPILEOPTS) $(F90COMPILEPATHS) $(F90COMPILECPPFLAGS) $<
       
      
# -----------------------------------------------------------------------------
shallow.exe: module_para.o module_array.o cs.o dif.o euler.o haurwitz.o  main.o
	$(F90LINKER) $(F90LINKOPTS) $(F90LINKPATHS) -o $@ $^ $(F90LINKLIBS)

# -----------------------------------------------------------------------------
.PHONY: dust clean distclean
dust:
	rm -f *.log 
clean:
	rm -f *.exe *.o *.mod
distclean: dust clean

run:
	make shallow.exe
	./shallow.exe 
