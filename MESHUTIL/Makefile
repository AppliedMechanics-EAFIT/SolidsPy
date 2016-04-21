# ============================================================================
# Name        : Makefile
# ============================================================================

.PHONY: all clean

.SUFFIXES: .for

FC=gfortran -g3 -ggdb
LIBS=-lpthread
FLAGS=
FLAGS+=$(INCS) $(LIBS)

SOURCES=mesher.for mtxutil.for meshutil.for

OBJECTS=$(SOURCES:.for=.o)

all: $(SOURCES) mesher

mesher: $(OBJECTS)
	$(FC) $(FLAGS) -o $@ $(OBJECTS)


.for.o:
	$(FC) $(INCS) -c $< -o $@

clean:
	rm -f mesher *.mod *.o

echo:
	echo $(OBJS)
