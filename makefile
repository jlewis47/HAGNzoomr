f90 = ifort
fpath=./utils_zoom_f90
vpath %.f90 $(fpath)



all: centre_grafic extract_grafic
	
centre_grafic: center_grafic_notinteractive.f90
	$(f90) $(fpath)/center_grafic_notinteractive.f90 -o ./center_grafic_notinteractive

extract_grafic: extract_grafic_notinteractive.f90
	$(f90) $(fpath)/extract_grafic_notinteractive.f90 -o ./extract_grafic_notinteractive

clean:
	rm -f extract_grafic_notinteractive center_grafic_notinteractive