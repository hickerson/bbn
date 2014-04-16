FABLEFILES=new123.f newint.f nuccom.f nucrat.f
FABLEFLAGS=--separate-cmn-hpp --namespace=bbn --no-fem-do-safe
FABLEOPTS=--common-equivalence-simple
FLAGS=-Ifable_sources/fable -I/usr/include/gsl -g -std=c++11
LIBS=-lgsl -lgslcblas -lm
SOURCES=bbn.cpp


all:

f2c:
	f2c -C++ -I2 -A -s -r -r8 *.f

fable:
	. fable_build/setpaths.sh
	fable.cout $(FABLEFILES) $(FABLEFLAGS) > bbn-fable.cpp
	./purify_fable.sh bbn-fable.cpp cmn.hpp

clean:
	rm -rf *.o bbn *.dat

install-fable:
	wget http://cci.lbl.gov/fable_bundles/current/fable_bundle.selfx
	perl fable_bundle.selfx
	. ./fable_build/

bbn:
	rm -f bbn
	g++ $(FLAGS) bbn.cpp $(LIBS) -o bbn



run:
	rm -f bbn
	g++ -g -o bbn -I$(CLIBS) bbn.cpp
	rm -f *.dat
	./bbn

bbn-fable:
	rm -f bbn
	g++ -o bbn -I$(CLIBS) bbn-fable.cpp

comments:
	--top-procedure=IDENTIFIER
	--include-guard-suffix=STRING
	--dynamic-parameter="int array_size=100"
	--fortran-file-comments
	--no-fem-do-safe      
	--arr-nd-size-max=INTEGER (default: 256)
	--inline-all          
	--common-equivalence-simple=STRING
	--namespace=NAMESPACE
