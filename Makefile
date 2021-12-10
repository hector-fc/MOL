EXEC = molpack
LDFLAGS =
G77FLAGS = -c  -O3 -ffast-math -Wall
OBJECTS = molpack.o
$(EXEC): $(OBJECTS)
	gfortran -o $@ $(OBJECTS) -lm

molpack.o: molpack.f
	gfortran $(G77FLAGS) $<

clean:
	rm -f *.o *.txt *.eps molpack molpack.out
