# Grab the targets and sources as two batches
SOURCES = $(wildcard src/*.cxx)
OBJECTS = $(patsubst src%.cxx,build%.o,$(SOURCES))

# Set the directory of the FID headers and lib
FID_DIR = ../cpp

FLAGS = -std=c++11 -g -I$(FID_DIR)/include
LIBS = -lfftw3 -lm -L/usr/local/lib/ -lboost_filesystem -lboost_system -L$(FID_DIR)/lib/ -lFid

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