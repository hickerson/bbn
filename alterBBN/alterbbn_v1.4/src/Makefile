
.KEEP_STATE:

.SUFFIXES:	.o .c .h
.PRECIOUS:	.c .h libbbn.a

include FlagsForMake

.c.o:
	$(CC) -c $(CFLAGS) $<
.c.a:
	$(CC) -c $(CFLAGS) $<
	$(AR) r $@ $*.o; rm $*.o

all: libbbn.a
	@case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac
	   
clean:
	rm -f *.a

distclean: 
	rm -f *.a *.o *.x
	
libbbn.a: libbbn.a(general.o) libbbn.a(omega.o) libbbn.a(bbn.o) libbbn.a(bbnrate.o)
	$(RANL)

libbbn.a(general.o): general.c include.h
libbbn.a(omega.o): omega.c include.h
libbbn.a(bbn.o): bbn.c include.h
libbbn.a(bbn.o): bbnrate.c include.h
