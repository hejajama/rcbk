CXXFLAGS = `gsl-config --cflags` -g -Wall -pedantic -I ../amplitudelib #-fopenmp# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm  

include filelist.m

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
