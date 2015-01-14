CC=gcc
CFLAGS=-I.
LFLAGS=-lf2c -lm
DEPS=
BINS=new123
OBJS=nuccom.o nucrat.o newint.o new123.o
DATS=newint.dat new123.dat

.PHONY: clean

all: clean $(BINS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

new123: $(OBJS)
	gcc -o new123 $(OBJS) $(LFLAGS)

clean:
	rm -f $(BINS) $(OBJS) $(DATS)
