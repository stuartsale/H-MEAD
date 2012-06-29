# Makefile for MEAD-BB
# Written by Stuart Sale 17/11/10
RUN_DIR=../../
flags= -g -fopenmp
linking_flags= -lgsl -lgslcblas -lm -lprofiler -fopenmp -g

MEAD: bin_obj.o iso_obj.o helper.o iphas_obj.o mead.o sl_obj.o cat_read.o
	g++ -o H-MEAD bin_obj.o iso_obj.o helper.o iphas_obj.o sl_obj.o cat_read.o mead.o $(linking_flags)
	cp H-MEAD $(RUN_DIR)

bin_obj.o: bin_obj.cpp bin_obj.h
	g++ -c $(flags) bin_obj.cpp

cat_read.o: cat_read.cpp cat_read.h iphas_obj.h
	g++ -c $(flags) cat_read.cpp	

iso_obj.o: iso_obj.cpp iso_obj.h
	g++ -c $(flags) iso_obj.cpp

helper.o: helper.cpp helper.h
	g++ -c $(flags) helper.cpp

iphas_obj.o: iphas_obj.cpp bin_obj.h iso_obj.h helper.h iphas_obj.h
	g++ -c $(flags) iphas_obj.cpp

sl_obj.o: sl_obj.cpp sl_obj.h helper.h
	g++ -c $(flags) sl_obj.cpp

mead.o: mead.cpp bin_obj.h iso_obj.h helper.h iphas_obj.h sl_obj.h
	g++ -c $(flags) mead.cpp



