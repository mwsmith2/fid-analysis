TARGETS = $(patsubst prog/%.cxx,%,$(wildcard prog/*.cxx))
SOURCES = $(wildcard src/*.cxx)
SOURCES += $(wildcard prog/*.cxx)
HEADERS = $(SOURCES:.cxx=.hh)
OBJECTS = $(patsubst src%.cxx,build%.o,$(SOURCES))

FLAGS = -std=c++11 -O3 -Iinclude
LIBS = -lfftw3 -lm

FLAGS += $(shell root-config --cflags)
LIBS  += $(shell root-config --libs)

CXX = clang++

all:

build/%.o: src/%.cxx
	$(CXX) $(FLAGS) -o $@ -c $< $(LIBS)

$(TARGETS): $(OBJECTS)
	$(CXX) $(FLAGS) -o $@ $+ $(LIBS)

clean:
	rm -f $(TARGETS) build/* 