#-------------------------------------------------------------------
#
# Fortran
#
FFILES = new123.f newint.f nuccom.f nucrat.f
FFLAGS = --separate-cmn-hpp --namespace=bbn --no-fem-do-safe
FOPTS = --common-equivalence-simple


#-------------------------------------------------------------------
#
# commands
#
CC = colorgcc
RM = rm -rf
FCC = fable.out


#-------------------------------------------------------------------
#
# libraries 
#
GSLLIBS = -lgsl -lgslcblas
LIBS = -lm $(GSLLIBS)


#-------------------------------------------------------------------
#
# sources
#
MAIN = bbn
SRCS = bbn.cpp EvolutionParamaters.cpp
INCS = BigBangModel.hpp EvolutionParamaters.hpp
OBJS = $(SRCS:.cpp=.o)


#-------------------------------------------------------------------
#
# flags
#
DEBUG = -g -Wall
CFLAGS = $(DEBUG) -std=c++11 -fpermissive
IFLAGS = -Ifable_sources/fable -I/usr/include/gsl
LFLAGS = 


#-------------------------------------------------------------------
#
# arguments
#

default: $(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(IFLAGS) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c $< -o $@

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

bbn-fable:
	rm -f bbn
	g++ -o bbn -I$(CLIBS) bbn-fable.cpp

.PHONY: clean
clean:
	$(RM) *.o $(MAIN) *.dat

#depend: $(SRCS)
#	makedepend $(IFLAGS) $^

# NOT DELETE THIS LINE -- make depend needs it

