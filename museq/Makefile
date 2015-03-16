###################################################
#
# Makefile for pybam.so
#
###################################################

#
# Macros
#

BOOST_ROOT ?= share/apps/build/boost_1_51_0

INCLUDE = -I. 
CC = g++
CXX = g++
CXXFLAGS = -O3 -ggdb -DDEBUG_CHECKS $(INCLUDE)
LDFLAGS = -O3 -ggdb 
LDLIBS = 

BINS = \
		pybam.so

SOURCE = \
		src/pybam.cpp

#
# Build
#

all : $(BINS)

pybam.so : $(SOURCE)
		python setup.py --boost_source=$(BOOST_ROOT) install --prefix=./ --install-platlib=./

clean : 
		rm -f *.o *.so
		rm -rf build