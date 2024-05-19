F77 = mpif77
FFLAGS = -O3 -g #-std=legacy

olist = dpmeps.o enorm.o lmder1.o lmdif1.o lmpar.o qrsolv.o dpmpar.o \
	fdjac2.o lmder.o lmdif.o  qrfac.o src.o \
	hybrd1.o hybrd.o qform.o r1mpyq.o r1updt.o dogleg.o fdjac1.o \
	hybrj1.o hybrj.o

libminpack.a: $(olist)
	ar rv $@ $(olist)
	ranlib $@

.f.o:
	$(F77) $(FFLAGS) -c $*.f

dpmeps.o: toms-768/Src/Dp/dpmeps.f
	$(F77) $(FFLAGS) -c toms-768/Src/Dp/dpmeps.f

src.o: toms-768/Src/Dp/src.f
	$(F77) $(FFLAGS) -c toms-768/Src/Dp/src.f

clean:
	rm -f *.o libminpack.a

FILES = *.f *.f.orig Makefile README \
	toms-768/Drivers/Dp/*.f \
	toms-768/Drivers/Dp/RES? \
	toms-768/Src/Dp/*.f \
	toms-768/Doc/makefile

minpack.tar.bz2: $(FILES)
	@tar cjpvf $@ $^
