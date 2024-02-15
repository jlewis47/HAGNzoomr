f90 = ifort

path_code=./utils_zoom_f90

all:
	$(f90) $(path_code)/center_grafic_notinteractive.f90 -o ./center_grafic_notinteractive
	$(f90) $(path_code)/extract_grafic_notinteractive.f90 -o ./extract_grafic_notinteractive