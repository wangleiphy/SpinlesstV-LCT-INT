include $(ALPS_ROOT)/share/alps/include.mk
include compiler.flags

all: main

OBJ = io.o interaction_expansion.o fastupdate.o model.o solver.o observables.o measure.o main.o 

main: $(OBJ) 
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	$(RM) $(RMFLAGS) *.o

