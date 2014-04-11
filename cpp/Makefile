# Grab the targets and sources as two batches
TARGETS = $(patsubst prog/%.cxx,%,$(wildcard prog/*.cxx))
SOURCES = $(wildcard src/*.cxx)
HEADERS = $(SOURCES:.cxx=.hh)
OBJECTS = $(patsubst src%.cxx,build%.o,$(SOURCES))
# Grab the sources in prog
#SOURCES = $(wildcard prog/*.cxx)
#MODULES = $(patsubst prog%.cxx,build%.o,$(SOURCES))

FLAGS = -std=c++11 -O3 -Iinclude
LIBS = -lfftw3 -lm

FLAGS += $(shell root-config --cflags)
LIBS  += $(shell root-config --libs)
ARFLAGS = rvs 
LIBNAME = lib/libFid.a

CXX = clang++
AR = ar

all:

build/%.o: src/%.cxx
	$(CXX) $(FLAGS) -o $@ -c $<

lib: $(OBJECTS)
	$(AR) $(ARFLAGS) $(LIBNAME) $+

clean:
	rm -f $(TARGETS) build/* 