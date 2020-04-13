#
# Makefile for CHAMP
#
#           (C) A.Baldwin
#           University of Oxford
#           andrew.baldwin@chem.ox.ac.uk
#

CC = g++-8


LDFLAGS = -L/usr/local/lib -lgsl -lgslcblas  -lm -o 

VERSION = "7.0"

#UNAME_P = $(shell uname -p)
UNAME_S = $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	#CC = g++-8
        CC = g++
        MPICC = mpic++
	#CFLAGS :=  -DMAC_OS_X  -std=c++11 -O3 -Ofast -ffast-math -ftree-vectorize -Wall -Wno-switch-bool -static
	#CCFLAGS := -O3  -std=c++11 -DMAC_OS_X  -Wno-string-compare -Wno-return-type -Wno-writable-strings  -Wno-comment -static
	CCFLAGS := -O3 -std=c++11 -ffast-math -DMAC_OS_X  -Wno-string-compare -Wno-return-type -Wno-writable-strings  -Wno-comment -Wno-deprecated -Wno-unused-command-line-argument -Wno-dangling-else -static -I/usr/local/include -L/usr/local/lib

endif

ifeq ($(UNAME_S),Linux)
        CC = g++
        MPICC = /usr/lib64/openmpi/bin/mpic++
	CCFLAGS := -O3  -std=c++11 -DLINUX  -Wno-string-compare -Wno-return-type -Wno-writable-strings -static
	#CCFLAGS :=  -DLINUX -DNDEBUG -std=c++11 -O3 -Ofast -ffast-math -ftree-vectorize -Wno-comment  -Wno-write-strings  -static
endif


SOURCE1 = disty.c
#SOURCE2 = mcesCore_parallel.cpp
OBJECT1 = ${SOURCE1:.c=.o}
#OBJECT2 = ${SOURCE2:.cpp=.o}


clean:
	@echo "Removing objects" 
	rm -f $(OBJECT1) 


install:
	@echo "Compiling binaries"
	$(CC) -c $(SOURCE1) -o $(OBJECT1) $(CCFLAGS)
	$(CC) $(OBJECT1) $(LDFLAGS)  ../bin/champ_$(UNAME_S)

	#$(MPICC) -c $(SOURCE2) -o $(OBJECT2) $(CCFLAGS)
	#$(MPICC) $(OBJECT2) $(LDFLAGS)  ../MAGMA_$(VERSION)/bin/mcesCore_parallel_$(UNAME_S)


	@echo "================================================================="
	@echo "   CHAMP - version "$(VERSION) 
	@echo "   Compiled successfully "	
	@echo "                    (C) CHAMP"
	@echo "                        A.Baldwin"
	@echo "                        University of Oxford  "          
	@echo "                        andrew.baldwin@chem.ox.ac.uk "         
	@echo "================================================================="
