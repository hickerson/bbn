FORTRAN=new123.f nuccom.f newint.f nucrat.f
FABLEFLAGS=--namespace=bbn --no-fem-do-safe --separate-cmn-hpp

all: c

f2c:
	f2c -C++ -I2 -A -s -r -r8 *.f

fable:
	. fable_build/setpaths.sh
	fable.cout $(FORTRAN) $(FABLEFLAGS) > bbn.cpp
	. ./purify_fable.sh

noc:
	rm -rf *.cpp *.hpp

install-fable:
	wget http://cci.lbl.gov/fable_bundles/current/fable_bundle.selfx
	perl fable_bundle.selfx
	. ./fable_build/

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
