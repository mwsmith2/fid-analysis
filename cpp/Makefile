# Grab the targets and sources as two batches
SOURCES = $(wildcard src/*.cxx)
HEADERS = $(wildcard include/*.h)
OBJECTS = $(patsubst src%.cxx,build%.o,$(SOURCES))

# Versioning info
MAJOR=0
MINOR=6.1
SONAME=libfid.so
LIBNAME=$(SONAME).$(MAJOR).$(MINOR)
PREFIX=/usr/local

# Some optional flags
DEBUG = -g -pg
OPTIMIZE = -O3

# Figure out the architecture
UNAME_S = $(shell uname -s)

# Clang compiler
ifeq ($(UNAME_S), Darwin)
	CXX = clang++
	CC  = clang
	FLAGS = -std=c++11
	LDCONFIG = cp $(PREFIX)/lib/$(LIBNAME) $(PREFIX)/lib/$(SONAME).$(MAJOR)
endif

# Gnu compiler
ifeq ($(UNAME_S), Linux)
	CXX = g++
	CC  = gcc
	FLAGS = -std=c++0x 
	LDCONFIG = ldconfig -n -v $(PREFIX)/lib
endif

FLAGS += -Wall -fPIC $(DEBUG) $(OPTIMIZE) -Iinclude
LIBS = -lfftw3 -lm

FLAGS += $(shell root-config --cflags)
LIBS  += $(shell root-config --libs)

all: $(LIBNAME)

build/%.o: src/%.cxx
	$(CXX) $(FLAGS) -o $@ -c $<

$(LIBNAME): $(OBJECTS)
	$(CXX) -shared -fPIC $+ -o lib/$(LIBNAME) $(LIBS)

install:
	mkdir -p $(PREFIX)/lib
	cp lib/$(LIBNAME) $(PREFIX)/lib
	mkdir -p $(PREFIX)/include
	cp $(HEADERS) $(PREFIX)/include
	ln -sf $(PREFIX)/lib/$(LIBNAME) $(PREFIX)/lib/$(SONAME)
	$(LDCONFIG)

uninstall:
	rm -f $(PREFIX)/lib/$(SONAME)*
	rm -f $(patsubst include/%,$(PREFIX)/include/%,$(HEADERS))

clean:
	rm -f $(TARGETS) build/* 
