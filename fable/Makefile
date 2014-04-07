FABLEFILES=new123.f newint.f nuccom.f nucrat.f
FABLEFLAGS=--namespace=bbn --no-fem-do-safe --separate-cmn-hpp
FABLEOPTS=--common-equivalence-simple

all: c

f2c:
	f2c -C++ -I2 -A -s -r -r8 *.f

fable:
	. fable_build/setpaths.sh
	fable.cout $(FABLEFILES) $(FABLEFLAGS) > bbn.cpp
	. ./purify_fable.sh

noc:
	rm -rf *.cpp *.hpp

install-fable:
	wget http://cci.lbl.gov/fable_bundles/current/fable_bundle.selfx
	perl fable_bundle.selfx
	. ./fable_build/

bbn:
	g++ -o bbn -Ifable_sources/fable bbn.cpp

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
