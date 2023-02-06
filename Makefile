CXXFLAGS = `gsl-config --cflags` -O3 -openmp -pedantic 
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
	rm -f $(OBJECTS) $(AMPLITUDELIBO)
	rm -f rbk	
