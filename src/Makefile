# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include template


.DEFAULT:
	-touch $@
all: adr2d
adr2d.o: ./adr2d.f90 decimal.o tipos.o storage.o boundary.o system.o plot.o util.o error.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./adr2d.f90
boundary.o: ./boundary.f90 decimal.o funciones.o tipos.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./boundary.f90
c_fortran_dgssv_crs.o: ./c_fortran_dgssv_crs.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	./c_fortran_dgssv_crs.c
decimal.o: ./decimal.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./decimal.f90
error.o: ./error.f90 decimal.o tipos.o quadrature.o funciones.o util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./error.f90
funciones.o: ./funciones.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./funciones.f90
plot.o: ./plot.f90 decimal.o tipos.o util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./plot.f90
quadrature.o: ./quadrature.f90 tipos.o decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./quadrature.f90
storage.o: ./storage.f90 decimal.o tipos.o quadrature.o funciones.o util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./storage.f90
system.o: ./system.f90 decimal.o tipos.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./system.f90
tipos.o: ./tipos.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./tipos.f90
util.o: ./util.f90 decimal.o tipos.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./util.f90
SRC = ./boundary.f90  ./util.f90 ./tipos.f90 ./decimal.f90 ./c_fortran_dgssv_crs.c ./funciones.f90 ./adr2d.f90 ./system.f90 ./error.f90 ./quadrature.f90 ./plot.f90 ./storage.f90
OBJ = boundary.o util.o tipos.o decimal.o c_fortran_dgssv_crs.o funciones.o  adr2d.o system.o error.o quadrature.o plot.o storage.o
clean: neat
	-rm -f .cppdefs $(OBJ) ../bin/adr2d *.mod
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
adr2d: $(OBJ) 
	$(LD) $(OBJ) -o ../bin/adr2d  $(LDFLAGS)
