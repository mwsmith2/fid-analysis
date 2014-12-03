# Grab the targets and sources as two batches
SOURCES = $(wildcard src/*.cxx)
TARGETS = $(patsubst src/%.cxx,bin/%,$(SOURCES))

# Debugging/profiling flags
DEBUG = -g -pg

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
LIBS = -lfftw3 -lm -lboost_filesystem -lboost_system -lfid

FLAGS += $(shell root-config --cflags)
LIBS  += $(shell root-config --libs)

all: $(TARGETS)

build/%.o: src/%.cxx
	$(CXX) $(FLAGS) -o $@ -c $< 

bin/%: build/%.o
	$(CXX) $(FLAGS) -o $@ $+ $(LIBS)

clean:
	rm -f build/* 