CC = mpicc
CFLAGS = -g -Wall -MMD
LDFLAGS = -lm

.PHONY: clean all

all: raytracer

clean:
	rm *.bmp
	rm *.d
	rm *.o
	rm raytracer

raytracer: vector.o bmp.o triangle.o textProcessor/textProcessor.o

-include *.d
