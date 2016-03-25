.KEEP_STATE:

#
VERSION = v1.7
YEAR = 2016
AUTHORS = Dr. Kevin Peter Hickerson

# Choose your compilers here (usually gcc/gfortran on Linux systems):
CC = gcc
CXX = g++
CFLAGS=-Wmaybe-uninitialized -Wuninitialized -Winit-self -g -std=gnu++11 -fext-numeric-literals

# CC = clang
# CFLAGS= -O3

# CC = icc
# CFLAGS = -O3 -Wall

MAKE = make
AR = ar

.SUFFIXES:	.o .cpp .h
.PRECIOUS:	.cpp .h libbbn.a

CINCLUDE= -I./src -L./src

all: libbbn.a
	@case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac;\
	echo ' ';\
   	echo 'Please run "make name" to compile "name.cpp" into "name.x"';\
	echo ' '

%.cpp:: %.cpp libbbn.a
	$(CC) -c $(CFLAGS) $@;
	$(CXX) -o $*.x $(CFLAGS) $(CINCLUDE) $*.o -lbbn -lm;
	@rm -f $*.o;
	@touch $*.x

%:: %.cpp libbbn.a
	$(CC) -c $(CFLAGS) $*.cpp;
	$(CXX) -o $*.x $(CFLAGS) $(CINCLUDE) $*.o -lbbn -lm;
	@rm -f $*.o;
	@touch $*.x

clean:
	rm -rf tmp *.x;
	@echo > src/FlagsForMake;
	$(MAKE) -C src/ clean
	
distclean: 
	rm -rf tmp *.x;
	@echo > src/FlagsForMake;
	$(MAKE) -C src/ distclean
	
libbbn.a: 
	@echo;
	@echo AlterBBN v1.4 - A. Arbey 2013;
	@echo AlterBBN $(VERSION) - K. P. Hickerson $(YEAR);
	@echo;
	@echo CC = $(CC) > src/FlagsForMake;\
	echo CXX = $(CXX) > src/FlagsForMake;\
	echo CFLAGS = $(CFLAGS) >> src/FlagsForMake;\
	echo MAKE = $(MAKE) >> src/FlagsForMake;\
	echo AR = $(AR) >> src/FlagsForMake;
	$(MAKE) -C src/ libbbn.a

save: 
	rm -f alterbbn_$(VERSION).tgz;\
	mkdir alterbbn_$(VERSION);\
	cp -p README alterbbn_$(VERSION)/;\
	cp -p stand_cosmo.c alterbbn_$(VERSION)/;\
	cp -p alter_eta.c alterbbn_$(VERSION)/;\
	cp -p alter_etannutau.c alterbbn_$(VERSION)/;\
	cp -p alter_neutrinos.c alterbbn_$(VERSION)/;\
	cp -p alter_reheating.c alterbbn_$(VERSION)/;\
	cp -p alter_standmod.c alterbbn_$(VERSION)/;\
	cp -p alter_fierz.cpp alterbbn_$(VERSION)/;\
	cp -p Makefile alterbbn_$(VERSION)/;\
	mkdir alterbbn_$(VERSION)/src;\
	cp -p src/*.h alterbbn_$(VERSION)/src/;\
	cp -p src/*.cpp alterbbn_$(VERSION)/src/;\
	cp -rp src/sgStar_heff alterbbn_$(VERSION)/src/;\
	cp -p src/Makefile alterbbn_$(VERSION)/src/;\
	tar czvf alterbbn_$(VERSION).tgz alterbbn_$(VERSION);\
	rm -rf alterbbn_$(VERSION)
