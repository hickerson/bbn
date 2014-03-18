all: c

f2c:
	f2c -C++ -I2 -A -s -r -r8 *.f

fable:
	fable

noc:
	rm -rf *.c

install-fable:
	wget http://cci.lbl.gov/fable_bundles/current/fable_bundle.selfx
	perl fable_bundle.selfx
	./fable_install_script.csh
