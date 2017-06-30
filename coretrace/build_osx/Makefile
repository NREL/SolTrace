VPATH = ..
CC = gcc
CXX = g++
CFLAGS = -fPIC -Wall -g -O3 -I../ 
CXXFLAGS = $(CFLAGS)

OBJECTS = \
	treemesh.o \
	apertureplane.o \
	determineelementintersectionnew.o \
	dumpsys.o \
	errors.o \
	generateray.o \
	input.o \
	interaction.o \
	intersect.o \
	mathproc.o \
	newzstartforcubicsplinesurf.o \
	quadricsurfaceclosedform.o \
	root432.o \
	spencerandmurtysurfaceclosedform.o \
	stapi.o \
	suntoprimarystage.o \
	surface.o \
	surfacenormalerrors.o \
	surfacezatxypair.o \
	torusclosedform.o \
	raytrace.o \
	types.o \
	vshot.o


TARGET=coretrace.a

$(TARGET):$(OBJECTS)
	ar rs $(TARGET) $(OBJECTS)

clean:
	rm -rf $(TARGET) $(OBJECTS)
