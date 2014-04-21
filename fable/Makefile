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

#
# libraries 
#
LIBS = -lgsl -lgslcblas -lm

#
# sources
#
MAIN = bbn
SRCS = bbn.cpp EvlolutionParamaters.cpp
INCS = BigBangModel.hpp EvlolutionParamaters.hpp
OBJS = $(SRCS:.cpp=.o)

#
# options
#
DEBUG = -g -Wall
CFLAGS = $(DEBUG) -std=c++11 -fpermissive 
IFLAGS = -Ifable_sources/fable -I/usr/include/gsl
LFLAGS = 



all: $(MAIN)
	@echo "Make all..."

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(IFLAGS) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c $< -o $@

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


.PHONY: clean
clean:
	$(RM) *.o $(MAIN) *.dat

#depend: $(SRCS)
#	makedepend $(IFLAGS) $^

# NOT DELETE THIS LINE -- make depend needs it

