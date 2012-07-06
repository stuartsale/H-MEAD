# Makefile for MEAD-BB
# Written by Stuart Sale 17/11/10
HOST=$(shell hostname)
ifeq ($(HOST),orion)
	RUN_DIR=../../
	flags= -g -fopenmp
	linking_flags= -lgsl -lgslcblas -lm -lprofiler -fopenmp -g
	gsl_lib=/usr/lib/
	gsl_include=/usr/include/
else 
	ifeq ($(HOST),newhydra.physics.ox.ac.uk)
		RUN_DIR=../../
		flags= -g -fopenmp
		linking_flags= -lgsl -lgslcblas -lm -fopenmp -g
		gsl_include=/usr/local/shared/gsl-1.12/include/
		gsl_lib=/usr/local/shared/gsl-1.12/lib/
	else	#defaults
		RUN_DIR=../../
		flags= -g -fopenmp
		linking_flags= -lgsl -lgslcblas -lm -fopenmp -g
		gsl_lib=/usr/lib/
		gsl_include=/usr/include/
	endif

endif



MEAD: bin_obj.o iso_obj.o helper.o iphas_obj.o mead.o sl_obj.o cat_read.o
	g++ -o H-MEAD bin_obj.o iso_obj.o helper.o iphas_obj.o sl_obj.o cat_read.o mead.o -I$(gsl_include) -L$(gsl_lib) $(linking_flags)
	cp H-MEAD $(RUN_DIR)

bin_obj.o: bin_obj.cpp bin_obj.h
	g++ -c $(flags) bin_obj.cpp -I$(gsl_include) -L$(gsl_lib)

cat_read.o: cat_read.cpp cat_read.h iphas_obj.h
	g++ -c $(flags) cat_read.cpp -I$(gsl_include) -L$(gsl_lib)

iso_obj.o: iso_obj.cpp iso_obj.h
	g++ -c $(flags) iso_obj.cpp -I$(gsl_include) -L$(gsl_lib)

helper.o: helper.cpp helper.h
	g++ -c $(flags) helper.cpp -I$(gsl_include) -L$(gsl_lib)

iphas_obj.o: iphas_obj.cpp bin_obj.h iso_obj.h helper.h iphas_obj.h
	g++ -c $(flags) iphas_obj.cpp -I$(gsl_include) -L$(gsl_lib)

sl_obj.o: sl_obj.cpp sl_obj.h helper.h
	g++ -c $(flags) sl_obj.cpp  -I$(gsl_include) -L$(gsl_lib)

mead.o: mead.cpp bin_obj.h iso_obj.h helper.h iphas_obj.h sl_obj.h
	g++ -c $(flags) mead.cpp -I$(gsl_include) -L$(gsl_lib)



