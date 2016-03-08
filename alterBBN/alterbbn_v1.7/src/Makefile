
.KEEP_STATE:

.SUFFIXES:	.o .cpp .h
.PRECIOUS:	.cpp .h libbbn.a

include FlagsForMake

.cpp.o:
	$(CXX) -c $(CFLAGS) $<
.cpp.a:
	$(CXX) -c $(CFLAGS) $<
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

libbbn.a(general.o): general.cpp include.h
libbbn.a(omega.o): omega.cpp include.h
libbbn.a(bbn.o): bbn.cpp include.h
libbbn.a(bbn.o): bbnrate.cpp include.h
