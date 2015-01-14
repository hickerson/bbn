CC=gcc
CFLAGS=-I.
DEPS=

%.o: %.c $(DEPS)
    $(CC) -c -o $@ $< $(CFLAGS)

new123: nuccom.o nucrat.o newint.o new123.o
    gcc -o new123 nuccom.o nucrat.o newint.o new123.o -lf2c -lm

