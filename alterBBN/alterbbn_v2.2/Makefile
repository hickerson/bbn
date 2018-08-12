.KEEP_STATE:

#
VERSION = v2.0
YEAR = 2018

# Choose your compilers here (usually gcc/gfortran on Linux systems):
CC = gcc
CFLAGS= -O3 -pipe -fomit-frame-pointer
CFLAGS_MP= -fopenmp

#CC = clang
#CFLAGS= -O3 -pipe -fomit-frame-pointer -fopenmp
#CFLAGS_MP= -fopenmp

#CC = icc
#CFLAGS = -O3
#CFLAGS_MP= -qopenmp

MAKE = make
MAKE_MP = -j$(nproc)

AR = ar
ARFLAGS = rcsU
#ARFLAGS = rcs


.SUFFIXES:	.o .c .h
.PRECIOUS:	.c .h libbbn.a

CINCLUDE= -I./src -L./src

all: libbbn.a
	@case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) $(CFLAGS_MP) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac;\
	echo ' ';\
   	echo 'Please run "make name" to compile "name.c" into "name.x"';\
	echo ' '

%.c:: %.c libbbn.a
	$(CC) -c $(CFLAGS) $@;
	$(CC) -o $*.x $(CFLAGS) $(CFLAGS_MP) $(CINCLUDE) $*.o -lbbn -lm;
	@rm -f $*.o;
	@touch $*.x

%:: %.c libbbn.a
	$(CC) -c $(CFLAGS) $*.c;
	$(CC) -o $*.x $(CFLAGS) $(CFLAGS_MP) $(CINCLUDE) $*.o -lbbn -lm;
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
	@echo AlterBBN $(VERSION) - A. Arbey, K. Hickerson, E. Jenssen $(YEAR);
	@echo;
	@echo CC = $(CC) > src/FlagsForMake;\
	echo CFLAGS = $(CFLAGS) >> src/FlagsForMake;\
	echo CFLAGS_MP = $(CFLAGS_MP) >> src/FlagsForMake;\
	echo MAKE = $(MAKE) >> src/FlagsForMake;\
	echo AR = $(AR) >> src/FlagsForMake;\
	echo ARFLAGS = $(ARFLAGS) >> src/FlagsForMake;
	$(MAKE) $(MAKE_mp) -C src/ libbbn.a

save: 
	rm -f alterbbn_$(VERSION).tgz;\
	mkdir alterbbn_$(VERSION);\
	cp -p README alterbbn_$(VERSION)/;\
	cp -p stand_cosmo.c alterbbn_$(VERSION)/;\
	cp -p alter_eta.c alterbbn_$(VERSION)/;\
	cp -p alter_etannutau.c alterbbn_$(VERSION)/;\
	cp -p alter_neutrinos.c alterbbn_$(VERSION)/;\
	cp -p alter_phi.c alterbbn_$(VERSION)/;\
	cp -p alter_reheating.c alterbbn_$(VERSION)/;\
	cp -p alter_standmod.c alterbbn_$(VERSION)/;\
	cp -p alter_wimps.c alterbbn_$(VERSION)/;\
	cp -p Makefile alterbbn_$(VERSION)/;\
	mkdir alterbbn_$(VERSION)/src;\
	cp -p src/*.h alterbbn_$(VERSION)/src/;\
	cp -p src/*.c alterbbn_$(VERSION)/src/;\
	cp -rp src/sgStar_heff alterbbn_$(VERSION)/src/;\
	cp -p src/Makefile alterbbn_$(VERSION)/src/;\
	tar czvf alterbbn_$(VERSION).tgz alterbbn_$(VERSION);\
	rm -rf alterbbn_$(VERSION)
