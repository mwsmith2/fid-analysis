# Grab the targets and sources as two batches
SOURCES = $(wildcard src/*.cxx)
OBJECTS = $(patsubst src%.cxx,build%.o,$(SOURCES))
TARGETS = $(patsubst src/%.cxx,%,$(SOURCES))

# Debugging/profiling flags
DEBUG = -g -pg

# Set the directory of the FID headers and lib
FID_DIR = ../cpp

# Figure out the architecture
UNAME_S = $(shell uname -s)

# Clang compiler
ifeq ($(UNAME_S), Darwin)
	CXX = clang++
	CC  = clang
	FLAGS = $(DEBUG) #-std=c++11 $(DEBUG)
endif

# Gnu compiler
ifeq ($(UNAME_S), Linux)
	CXX = g++
	CC  = gcc
	FLAGS = -std=c++0x $(DEBUG)
endif

FLAGS += -I$(FID_DIR)/include
LIBS = -lfftw3 -lm -lboost_filesystem -lboost_system -L$(FID_DIR)/lib/ -lfid-shared

FLAGS += $(shell root-config --cflags)
LIBS  += $(shell root-config --libs)

all: $(TARGETS)

build/%.o: src/%.cxx
	$(CXX) $(FLAGS) -o $@ -c $< 

%: build/%.o
	$(CXX) $(FLAGS) -o build/$@ $+ $(LIBS)

clean:
	rm -f build/* 