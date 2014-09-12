# Grab the targets and sources as two batches
SOURCES = $(wildcard src/*.cxx)
HEADERS = $(SOURCES:.cxx=.hh)
OBJECTS = $(patsubst src%.cxx,build%.o,$(SOURCES))
# Grab the sources in prog
#SOURCES = $(wildcard prog/*.cxx)
#MODULES = $(patsubst prog%.cxx,build%.o,$(SOURCES))

# Some optional flags
DEBUG = -g -pg
#OPTIMIZE = -O3

# Figure out the architecture
UNAME_S = $(shell uname -s)

# Clang compiler
ifeq ($(UNAME_S), Darwin)
	CXX = clang++
	CC  = clang
	FLAGS = -std=c++11 $(DEBUG) $(OPTIMIZE)
endif

# Gnu compiler
ifeq ($(UNAME_S), Linux)
	CXX = g++
	CC  = gcc
	FLAGS = -std=c++0x $(DEBUG) $(OPTIMIZE)
endif

LIBNAME = lib/libFid.a
AR = ar
ARFLAGS = rvs 

FLAGS += -Iinclude
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