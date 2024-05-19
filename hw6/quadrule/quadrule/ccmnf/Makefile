# $Id: Makefile,v 1.4 2020/01/09 13:14:55 zlb Exp $

all: ccmnf ccmfr

CC = gcc -fPIC
CFLAGS = -Wall -O2 -g
LINKER = gcc

#############################################
ccmnf: test.o ccmNF.o
	$(LINKER) -o $@ $^ -lm

ccmNF.o: ccmNF.c ccmNF.h
	$(CC) $(CFLAGS) -c ccmNF.c

test.o: test.c
	$(CC) $(CFLAGS) -c test.c

#############################################

ccmfr: test2.o ccmNF.o
	$(LINKER) -o $@ $^ -lm

test2.o: test2.c
	$(CC) $(CFLAGS) -c test2.c

############################################

# clean object and executive files
clean:
	rm -rf *.o ccmnf ccmfr

FILES = ccmNF.c ccmNF.h Makefile test.c test2.c .cvsignore

ccmnf.zip: $(FILES)
	@zip -9 -u -y $@ $^
