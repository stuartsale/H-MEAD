# Makefile for MEAD-BB
# Written by Stuart Sale 17/11/10
HOST=$(shell hostname)
ifeq ($(HOST),orion)
	RUN_DIR=../../
	flags= -wd2196 -openmp -O3
	linking_flags= -lgsl -lgslcblas -lm -lprofiler -openmp -lCCfits -mkl
	gsl_lib=/usr/lib/
	gsl_include=/usr/include/
	CCfits_include=/usr/include/
	CCfits_lib=/usr/lib/x86_64-linux-gnu/
	EIGEN_include=/home/stuart/work/soft9/eigen/
	MKL_include=/opt/intel/composer_xe_2013.5.192/mkl/include/
	COMPILER=icpc
	OMPI_CXX :=icpc
else 
	ifeq ($(HOST),newhydra.physics.ox.ac.uk)
		RUN_DIR=../../
		flags= -wd2196 -openmp -O3
		linking_flags= -lgsl -lgslcblas -lm -openmp  -lCCfits -mkl
		gsl_include=/usr/local/shared/gsl-1.12/include/
		gsl_lib=/usr/local/shared/gsl-1.12/lib/
		CCfits_include=/usersVol1/sale/soft9/include/
		CCfits_lib=/usersVol1/sale/soft9/lib/
		EIGEN_include=/usersVol1/sale/soft9/eigen/
		MKL_include=/opt/intel/composer_xe_2013.5.192/mkl/include/
		COMPILER=icpc
	else
		ifeq ($(HOST),hydra.physics.ox.ac.uk)
			RUN_DIR=../../
			flags= -openmp -O3
			linking_flags= -lgsl -lgslcblas -lm -openmp -lCCfits -mkl
			gsl_include=/usr/local/shared/gsl-1.12/include/
			gsl_lib=/usr/local/shared/gsl-1.12/lib/
			CCfits_include=/usersVol1/sale/soft9/include/
			CCfits_lib=/usersVol1/sale/soft9/lib/
			EIGEN_include=/usersVol1/sale/soft9/eigen/
			MKL_include=/opt/intel/composer_xe_2013.5.192/mkl/include/
			COMPILER=icpc
		else
			ifeq ($(HOST),auriga)
				RUN_DIR=../../
				flags= -wd2196 -openmp -O3
				linking_flags= -lgsl -lgslcblas -lm -openmp -lCCfits -mkl
				gsl_lib=/usr/lib/
				gsl_include=/usr/include/
				CCfits_include=/usr/include/
				CCfits_lib=/usr/lib/x86_64-linux-gnu/
				MKL_include=/opt/intel/composer_xe_2013.5.192/mkl/include/
				EIGEN_include=/home/stuart/work/soft9/eigen/
				COMPILER=icpc
				OMPI_CXX :=icpc
			else	#defaults
				RUN_DIR=../../
				flags= -g -fopenmp
				linking_flags= -lgsl -lgslcblas -lm -fopenmp -lCCfits -g -mkl
				gsl_lib=/usr/lib/
				gsl_include=/usr/include/
				CCfits_include=/usr/include/
				CCfits_lib=/usr/lib/x86_64-linux-gnu/
				MKL_include=/opt/intel/composer_xe_2013.5.192/mkl/include/
				COMPILER=icpc
			endif
		endif
	endif

endif





MEAD: bin_obj.o iso_obj.o helper.o iphas_obj.o mead.o sl_obj.o cat_read.o LF.o SFD_read.o
	$(COMPILER) -o H-MEAD bin_obj.o iso_obj.o helper.o iphas_obj.o sl_obj.o cat_read.o LF.o mead.o SFD_read.o -I$(gsl_include) -L$(gsl_lib) -I$(CCfits_include) -L$(CCfits_lib) $(linking_flags) 
	cp H-MEAD $(RUN_DIR)

bin_obj.o: bin_obj.cpp bin_obj.h
	$(COMPILER) -c $(flags) bin_obj.cpp -I$(gsl_include) -L$(gsl_lib) -I$(CCfits_include) -L$(CCfits_lib)

cat_read.o: cat_read.cpp cat_read.h iphas_obj.h
	$(COMPILER) -c $(flags) cat_read.cpp -I$(gsl_include) -L$(gsl_lib) -I$(CCfits_include) -L$(CCfits_lib)

iso_obj.o: iso_obj.cpp iso_obj.h
	$(COMPILER) -c $(flags) iso_obj.cpp -I$(gsl_include) -L$(gsl_lib) -I$(CCfits_include) -L$(CCfits_lib)

helper.o: helper.cpp helper.h bin_obj.h iso_obj.h SFD_read.h
	$(COMPILER) -c $(flags) helper.cpp -I$(gsl_include) -L$(gsl_lib) -I$(CCfits_include) -L$(CCfits_lib) -I$(MKL_include)

iphas_obj.o: iphas_obj.cpp bin_obj.h iso_obj.h helper.h iphas_obj.h SFD_read.h
	$(COMPILER) -c $(flags) iphas_obj.cpp -I$(gsl_include) -L$(gsl_lib) -I$(CCfits_include) -L$(CCfits_lib)

LF.o: LF.cpp LF.h helper.h iphas_obj.h bin_obj.h iso_obj.h SFD_read.h
	$(COMPILER) -c $(flags) LF.cpp -I$(gsl_include) -L$(gsl_lib) -I$(CCfits_include) -L$(CCfits_lib)

sl_obj.o: sl_obj.cpp sl_obj.h helper.h iphas_obj.h cat_read.h bin_obj.h iso_obj.h LF.h SFD_read.h
	$(COMPILER) -c $(flags) sl_obj.cpp  -I$(gsl_include) -L$(gsl_lib) -I$(CCfits_include) -L$(CCfits_lib)

mead.o: mead.cpp bin_obj.h iso_obj.h helper.h iphas_obj.h sl_obj.h SFD_read.h
	$(COMPILER) -c $(flags) mead.cpp -I$(gsl_include) -L$(gsl_lib) -I$(CCfits_include) -L$(CCfits_lib)

SFD_read.o: SFD_read.cpp SFD_read.h
	$(COMPILER) -c $(flags) SFD_read.cpp -I$(gsl_include) -L$(gsl_lib) -I$(CCfits_include) -L$(CCfits_lib)



