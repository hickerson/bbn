CC=gcc
CFLAGS=-I.
LFLAGS=-lf2c -lm
DEPS=
BINS=new123
OBJS=nuccom.o nucrat.o newint.o
DATS=newint.dat new123.dat bbnint.dat bbn123.dat

.PHONY: clean

all: clean $(BINS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

new123: $(OBJS) new123.o
	gcc -o new123 $(OBJS) new123.o $(LFLAGS)

clean:
	rm -f $(BINS) $(OBJS) $(DATS)
