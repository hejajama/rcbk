CXXFLAGS = `gsl-config --cflags` -g -Wall -pedantic #-fopenmp# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm  

SOURCES = src/main.cpp src/tools.cpp src/amplitude.cpp src/solver.cpp src/interpolation.cpp 

OBJECTS=$(SOURCES:.cpp=.o)

all: rbk

rbk: $(OBJECTS) 
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o rbk 
.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f bk	
