CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6dev`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6dev`

MIST1_3D: MIST1a.o
	$(CC) -o MIST1a MIST1a.o $(LDFLAGS)
MIST1_3D.o: MIST1a.cpp

clean:
	$(RM) *.o *~ MIST1a 

