#
# Makefile for Coag_aerosol 1.0
#
# Compiler: GNU g++
# 
# Author: Zhichao Li @ AAQRL, Washington University in St. Louis
#

# main program name
MPROG = coag

# headers
IDIR = ../include
# object files *.o
ODIR = ./obj
# system dependent (boost library path)
LDIR = "/mnt/p/boost/boost_1_63_0/"

# project dependent
_DEPS = consts.h coag_odes.h beta_df.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
_OBJS = consts.o coag_odes.o beta_df.o coag.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

CC = g++ -std=c++11
# header flag
CFLAGS = -I$(IDIR)
# library flag
LIBS = -I$(LDIR)

# compiling
$(ODIR)/%.o : %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

# linking
$(MPROG) : $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ $(MPROG)
