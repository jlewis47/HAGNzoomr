fpath=./utils_zoom_f90

# Ifort configuration
F90=ifort
CFLAGS=-Ofast -cpp

# Gfortran configuration
# F90=gfortran
# CFLAGS=-ffree-line-length-none -std=f2008 -Ofast -Wall -x f95-cpp-input -J..

vpath %.f90 $(fpath)

LFLAGS= #-g -check all -traceback -O0

MODOBJ=


all: centre_grafic extract_grafic extract_grafic_cuboid get_music_refmask clean_mod 
	
centre_grafic: center_grafic.f90
	$(F90) $(LFLAGS) $^ -o $@

extract_grafic: extract_grafic.f90
	$(F90) $(LFLAGS) $^ -o $@

extract_grafic_cuboid: extract_grafic_cuboid.f90
	$(F90) $(LFLAGS) $^ -o $@

get_music_refmask: $(fpath)/utils.o $(fpath)/get_music_refmask.f90
	$(F90) $(LFLAGS) $^ -o $@

%.o: %.f90
	$(F90) $(CFLAGS) -c $^ -o $@


clean_mod:
	rm *.mod

clean:
	rm -f extract_grafic extract_grafic_cuboid center_grafic get_music_refmask $(fpath)/*.o $(fpath)/*.mod