# ============================================================================
# Name        : Makefile
# ============================================================================

.PHONY: all clean

.SUFFIXES: .for

FC=gfortran -g3 -ggdb
LIBS=-lpthread
FLAGS=
FLAGS+=$(INCS) $(LIBS)

SOURCES=contour.for

OBJECTS=$(SOURCES:.for=.o)

all: $(SOURCES) contour

contour: $(OBJECTS)
	$(FC) $(FLAGS) -o $@ $(OBJECTS)


.for.o:
	$(FC) $(INCS) -c $< -o $@

clean:
	rm -f contour *.mod *.o

echo:
	echo $(OBJS)
