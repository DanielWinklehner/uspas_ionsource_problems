CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6dev`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6dev`

ecr: ecr.o
	$(CC) -o ecr ecr.o $(LDFLAGS)
ecr.o: ecr.cpp

analysis: analysis.o
	$(CC) -o analysis analysis.o $(LDFLAGS)
analysis.o: analysis.cpp

clean:
	$(RM) *.o analysis ecr
