include $(ALPS_ROOT)/share/alps/include.mk
include compiler.flags

all: worm main

OBJ = io.o interaction_expansion.o splines.o fastupdate.o model.o solver.o utils.o worm.o wormupdate.o main.o 

worm: $(OBJ) wormobservables.o 
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

main: $(OBJ) observables.o measure.o unequaltime.o
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	$(RM) $(RMFLAGS) *.o

