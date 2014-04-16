# Grab the targets and sources as two batches
SOURCES = $(wildcard src/*.cxx)
HEADERS = $(SOURCES:.cxx=.hh)
OBJECTS = $(patsubst src%.cxx,build%.o,$(SOURCES))
# Grab the sources in prog
#SOURCES = $(wildcard prog/*.cxx)
#MODULES = $(patsubst prog%.cxx,build%.o,$(SOURCES))

# Figure out the architecture
UNAME_S = $(shell uname -s)

ifeq ($UNAME_S, Darwin)
	# Clang compiler
	CXX = clang++
	CC  = clang
	FLAGS = -std=c++11
endif

ifeq ($UNAME_S, Linux)
	# Gnu compiler
	CXX = g++
	CC  = gcc
	FLAGS = -std=c++0x

LIBNAME = lib/libFid.a
AR = ar
ARFLAGS = rvs 

FLAGS = -O3 -Iinclude
LIBS = -lfftw3 -lm

FLAGS += $(shell root-config --cflags)
LIBS  += $(shell root-config --libs)

all:

build/%.o: src/%.cxx
	$(CXX) $(FLAGS) -o $@ -c $<

lib: $(OBJECTS)
	$(AR) $(ARFLAGS) $(LIBNAME) $+

clean:
	rm -f $(TARGETS) build/* 