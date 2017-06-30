VPATH = ../src

WEXDIR = ../../../wex
WEXLIB = $(WEXDIR)/wexosxwx3.a
LKDIR = ../../../lk
LKLIB = $(LKDIR)/lkosxwx3.a
CORETRACEDIR = ../../coretrace
CORETRACELIB = $(CORETRACEDIR)/build_osx/coretrace.a

CC = gcc
CXX = g++
CFLAGS = -g -I. -I$(CORETRACEDIR) -I$(WEXDIR)/include -I$(LKDIR)/include -I../tcs -DLK_USE_WXWIDGETS `wx-config-3 --cflags`
CXXFLAGS = $(CFLAGS)
LDFLAGS = $(WEXLIB) $(LKLIB) $(CORETRACELIB) `wx-config-3 --libs` `wx-config-3 --libs aui` `wx-config-3 --libs stc` `wx-config-3 --libs gl` -lm -lcurl -ldl

OBJECTS  = soltrace.o \
	elementlist.o \
	fluxmap.o \
	intersections.o \
	geometry.o \
	optics.o \
	raydata.o \
	script.o \
	sunshape.o \
	project.o \
	trace.o
	
TARGET = SolTrace

$(TARGET): $(OBJECTS) $(WEXLIB) 
	$(CXX) -g -o $@ $^ $(LDFLAGS)
	SetFile -a t $@
	mkdir -p $(TARGET).app
	mkdir -p $(TARGET).app/Contents
	mkdir -p $(TARGET).app/Contents/MacOS
	mkdir -p $(TARGET).app/Contents/Resources
	cp $(TARGET) ./$(TARGET).app/Contents/MacOS
	cp $(TARGET).icns ./$(TARGET).app/Contents/Resources
	cp Info-$(TARGET).plist ./$(TARGET).app/Contents/Info.plist
	echo "-n APPL????" >> ./$(TARGET).app/Contents/PkgInfo
	
clean:
	rm -f $(TARGET) $(OBJECTS)
	rm -rf $(TARGET).app
