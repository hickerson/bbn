#
# Fortran
#
FFILES = new123.f newint.f nuccom.f nucrat.f
FFLAGS = --separate-cmn-hpp --namespace=bbn --no-fem-do-safe
FOPTS = --common-equivalence-simple

#
# commands
#
CC = g++
RM = rm -rf
FCC = fable.out


# path to any includes 
INCLUDES = -Ifable_sources/fable -I/usr/include/gsl

#
# libraries 
#
LIBS = -lgsl -lgslcblas -lm

#
# sources
#
SRCS = bbn.cpp EvlolutionParamaters.cpp

#
# options
#
DEBUG = -g -Wall
CFLAGS = $(DEBUG) -std=c++11 -fpermissive 
#LFLAGS = 

TARGET = bbn
MAIN = $(TARGET)

# 
# Generic Make file
#

.PHONY: depend clean

all: $(MAIN)
	@echo "Make all..."

$(MAIN): $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

#f2c:
#	f2c -C++ -I2 -A -s -r -r8 *.f

fable:
	. fable_build/setpaths.sh
	$(FCC) $(FFILES) $(FFLAGS) > bbn-fable.cpp
	./purify_fable.sh bbn-fable.cpp cmn.hpp

install-fable:
	wget http://cci.lbl.gov/fable_bundles/current/fable_bundle.selfx
	perl fable_bundle.selfx
	. ./fable_build/

#bbn: $(OBJS)
#	rm -f bbn
#	$(CC) $(LFLAGS) $(OBJS) $(LIBS) -o bbn

run:
	rm -f bbn
	g++ -g -o bbn -I$(CLIBS) bbn.cpp
	rm -f *.dat
	./bbn

bbn-fable:
	rm -f bbn
	g++ -o bbn -I$(CLIBS) bbn-fable.cpp

clean:
	$(RM) *.o $(MAIN) *.dat

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# NOT DELETE THIS LINE -- make depend needs it

