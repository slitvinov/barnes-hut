.POSIX:
.SUFFIXES: .o .c
.PHONY: clean lint
include ../conf.mk

O = \
main.o\

L = libbarnes-hut.a
H = barnes-hut.h

$L: $O; ar rv $@ $O && ranlib $@

.c.o:; $(CC) -c $(CFLAGS) $< -o $@
include dep.mk

lint:; make CFLAGS='-Wall -Wextra -g -O2'
dep:; for i in *.c; do $(CPP) -MM $$i; done > dep.mk

clean:; rm -f -- $O $L
