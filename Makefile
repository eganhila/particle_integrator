CXX=g++
LFlAGS=-lhdf5
LIBS=-lhdf5

SRCS=main.cpp integrator.cpp sim_dat.cpp
OBJS=$(subst .cc,.o,$(SRCS))

all: pinter

pinter: $(OBJS)
	$(CXX) $(LFLAGS) -o pinter $(OBJS) $(LIBS) 

main.o: main.cpp sim_dat.h integrator.h

integrator.o: integrator.cpp integrator.h

sim_dat.o: sim_dat.h sim_dat.cpp

clean:
	$(RM) $(OBJS)
