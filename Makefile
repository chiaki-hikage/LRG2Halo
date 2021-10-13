F90	= ifort
FFLAGS 	= -O2 -cm -w -vec_report0 -sox -openmp -openmp_report0

DIR=/usr/local/fftw
FFTWDIR3=$(DIR)/fftw-3.2.2
FFTWDIR2=$(DIR)/fftw-2.1.5
HFLAGS = $(DIR)/Healpix_2.01/src/f90/mod
LIBCFITS= $(DIR)/cfitsio/libcfitsio.a
LIBHPIX= $(DIR)/Healpix_2.01/lib/libhealpix.a

CC	= mpic++
CFLAGS 	= -O2 -openmp -openmp_report0
INCDIRS_FFTW	= -I$(FFTWDIR3)/dft -g -I$(FFTWDIR3)/rdft  -g -I$(FFTWDIR3)/kernel -g -I$(FFTWDIR3) -g -I$(FFTWDIR3)/api -g -I/usr/local/include -I$(FFTWDIR3)
LIBS = -L/usr/local/lib -lrdft -ldft -lkernel -lapi -lm -lstdc++ -lfftw3

almtool=$(HFLAGS)/pix_tools.o

PROG	= main

.SUFFIXES: .cpp

.cpp.o	:
	$(CC) $(INCDIRS_FFTW) -c $<

default	: $(PROG)

cpara.o: cpara.f90
	$(F90) $(FFLAGS) -c $<

main	: main.o field.o ran2.o utils.o
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

main.o : main.h field.h utils.h fourier_ver3.h cpara.h

field.o : utils.h fourier_ver3.h cpara.h

test	:
	./main

clean:
	-rm -f $(PROG) *.o core
tidy:
	-rm -f $(PROG) 


