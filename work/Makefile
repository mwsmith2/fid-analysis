# Grab the targets and sources as two batches
SOURCES = $(wildcard src/*.cxx)
OBJECTS = $(patsubst src%.cxx,build%.o,$(SOURCES))

# Set the directory of the FID headers and lib
FID_DIR = ../cpp

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

FLAGS += -O2 -I$(FID_DIR)/include
LIBS = -lfftw3 -lm -lboost_filesystem -lboost_system -L$(FID_DIR)/lib/ -lFid

FLAGS += $(shell root-config --cflags)
LIBS  += $(shell root-config --libs)

CXX = clang++

all:

build/%.o: src/%.cxx
	$(CXX) $(FLAGS) -o $@ -c $< 

%: build/%.o
	$(CXX) $(FLAGS) -o build/$@ $+ $(LIBS)

clean:
	rm -f build/* 