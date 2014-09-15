# Grab the targets and sources as two batches
SOURCES = $(wildcard src/*.cxx)
HEADERS = $(SOURCES:.cxx=.hh)
OBJECTS = $(patsubst src%.cxx,build%.o,$(SOURCES))

# Versioning info
MAJOR=0
MINOR=3
SONAME=libfid.so
LIBNAME=$(SONAME).$(MAJOR).$(MINOR)
LIBPREFIX=/usr/local/lib

# Some optional flags
DEBUG = -g -pg
OPTIMIZE = -O3

# Figure out the architecture
UNAME_S = $(shell uname -s)

# Clang compiler
ifeq ($(UNAME_S), Darwin)
	CXX = clang++
	CC  = clang
	FLAGS = -std=c++11 $(DEBUG) $(OPTIMIZE)
	LDCONFIG = cp $(LIBPREFIX)/$(LIBNAME) $(LIBPREFIX)/$(SONAME).$(MAJOR)
endif

# Gnu compiler
ifeq ($(UNAME_S), Linux)
	CXX = g++
	CC  = gcc
	FLAGS = -std=c++0x $(DEBUG) $(OPTIMIZE)
	LDCONFIG = ldconfig -n -v $(LIBPREFIX)/$(LIBNAME)
endif

FLAGS += -Iinclude
LIBS = -lfftw3 -lm

FLAGS += $(shell root-config --cflags)
LIBS  += $(shell root-config --libs)

all: $(LIBNAME)

build/%.o: src/%.cxx
	$(CXX) $(FLAGS) -o $@ -c $<

$(LIBNAME): $(OBJECTS)
	$(CXX) -shared -fPIC $+ -o $(LIBNAME) $(LIBS)

install:
	cp $(LIBNAME) $(LIBPREFIX)/$(LIBNAME)
	$(LDCONFIG)
	ln -s $(LIBPREFIX)/$(LIBNAME) $(LIBPREFIX)/$(SONAME)

uninstall:
	rm -f $(LIBPREFIX)/$(SONAME)*

clean:
	rm -f $(TARGETS) build/* 