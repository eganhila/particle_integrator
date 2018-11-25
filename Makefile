CXX=g++
LDFlAGS=-std=c++11 #-lhdf5
LIBS=-std=c++11 -lhdf5 

src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)


pinter: $(obj)
	$(CXX)  $(LIBS) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	rm *.o pinter 
