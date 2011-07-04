# Makefile for MEAD-BB
# Written by Stuart Sale 17/11/10
RUN_DIR=../../
flags= -O3 -fopenmp
linking_flags= -lgsl -lgslcblas -lm -lprofiler -fopenmp

MEAD: bin_obj.o iso_obj.o helper.o iphas_obj.o mead.o
	g++ -o H-MEAD bin_obj.o iso_obj.o helper.o iphas_obj.o mead.o newran03/newran1.o newran03/newran2.o newran03/myexcept.o newran03/simpstr.o newran03/extreal.o $(linking_flags)
	cp H-MEAD $(RUN_DIR)

bin_obj.o: bin_obj.cpp bin_obj.h
	g++ -c $(flags) bin_obj.cpp	

iso_obj.o: iso_obj.cpp iso_obj.h
	g++ -c $(flags) iso_obj.cpp

helper.o: helper.cpp helper.h
	g++ -c $(flags) helper.cpp

iphas_obj.o: iphas_obj.cpp bin_obj.h iso_obj.h helper.h iphas_obj.h
	g++ -c $(flags) iphas_obj.cpp

mead.o: mead.cpp bin_obj.h iso_obj.h helper.h iphas_obj.h
	g++ -c $(flags) mead.cpp



