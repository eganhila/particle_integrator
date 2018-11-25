CXX=g++
CXXFLAGS=-std=c++11 #-lhdf5
LIBS=-std=c++11 -lhdf5 -lgtest -lpthread 
ADD_FLAGS=

src = $(wildcard *.cc)
obj = $(src:.cc=.o)


pinter: $(obj)
	$(CXX)  $(CXXFLAGS) -o $@ $^ $(LIBS)

ptest: ADD_FLAGS += -DGTEST
ptest: $(obj)
	$(CXX)  $(CXXFLAGS) $(ADD_FLAGS) -o $@ $^ $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(ADD_FLAGS) -o $@ -c $<


.PHONY: clean
clean:
	rm *.o pinter ptest
