VPATH = ../src/


WEXDIR = ../../../wex
WEXLIB = $(WEXDIR)/wexuxwx3.a
LKDIR = ../../../lk
LKLIB = $(LKDIR)/lkuxwx3.a
CORETRACEDIR = ../../coretrace
CORETRACELIB = $(CORETRACEDIR)/build_linux/coretrace.a

CC = gcc
CXX = g++

CFLAGS = -g -O2 -std=c++11 -I. -I$(LKDIR)/include -I$(WEXDIR)/include -I$(CORETRACEDIR) -DLK_USE_WXWIDGETS `wx-config-3 --cflags`

CXXFLAGS = $(CFLAGS)
LDFLAGS =  $(WEXLIB) $(LKLIB) $(CORETRACELIB) `wx-config-3 --libs stc` `wx-config-3 --libs aui` `wx-config-3 --libs gl` `wx-config-3 --libs` -lm -lcurl  -lfontconfig -ldl 


OBJECTS = \
	console.o\
	elementlist.o \
	fluxmap.o \
	intersections.o \
	geometry.o \
	optics.o \
	raydata.o \
	script.o \
	sunshape.o \
	soltrace.o \
	project.o \
	trace.o

TARGET = soltrace

$(TARGET): $(OBJECTS) $(WEXLIB)
	$(CXX) -g -o $@ $^ $(LDFLAGS)
	mkdir -p ../deploy/linux_64
	cp $(TARGET) ../deploy/linux_64

clean:
	rm -f $(TARGET) $(OBJECTS)

