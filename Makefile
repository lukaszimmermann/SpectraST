###################################################################################
# THE SOFTWARE IS PROVIDED BY THE INSTITUTE FOR SYSTEMS BIOLOGY (ISB)             #
# "AS IS" AND "WITH ALL FAULTS." ISB MAKES NO REPRESENTATIONS OR WARRANTI         #
# ES OF ANY KIND CONCERNING THE QUALITY, SAFETY OR SUITABILITY OF THE             #
# SOFTWARE, EITHER EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IM-     #
# PLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,          #
# OR NON-INFRINGEMENT.                                                            #
#                                                                                 #
# ISB MAKES NO REPRESENTATIONS OR WARRANTIES AS TO THE TRUTH, ACCURACY            #
# OR COMPLETENESS OF ANY STATEMENTS, INFORMATION OR MATERIALS CONCERNING          #
#                                                                                 #
# THE SOFTWARE THAT IS CONTAINED IN ANY DOCUMENTATION INCLUDED WITH THE           #
# SOFTWARE OR ON AND WITHIN ANY OF THE WEBSITES OWNED AND OPERATED BY ISB         #
#                                                                                 #
# IN NO EVENT WILL ISB BE LIABLE FOR ANY INDIRECT, PUNITIVE, SPECIAL,             #
# INCIDENTAL OR CONSEQUENTIAL DAMAGES HOWEVER THEY MAY ARISE AND EVEN IF          #
# ISB HAVE BEEN PREVIOUSLY ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.            #
#                                                                                 #
###################################################################################

SRC= src
SRCEXTERN=$(SRC)/extern
BUILD= build
BUILDEXTERN= $(BUILD)/extern
OBJ= obj

ARCHIVE=tar.gz
ARCHIVE_BOOST=tar.bz2

# DEPENDENCIES
# EXPAT
DEP_EXPAT_V=$(shell echo $(SRCEXTERN)/expat*.$(ARCHIVE) | grep -Eo '[0-9]+(.[0-9]+)*')
DEP_EXPAT_BUILD:=$(abspath $(BUILDEXTERN)/expat-$(DEP_EXPAT_V))
DEP_EXPAT_SRC=$(SRCEXTERN)/expat-$(DEP_EXPAT_V)

# zlib
DEP_ZLIB_V=$(shell echo $(SRCEXTERN)/zlib*.$(ARCHIVE) | grep -Eo '[0-9]+(.[0-9]+)*')
DEP_ZLIB_BUILD := $(abspath $(BUILDEXTERN)/zlib-$(DEP_ZLIB_V))
DEP_ZLIB_SRC := $(SRCEXTERN)/zlib-$(DEP_ZLIB_V)

# gsl
DEP_GSL_V=$(shell echo $(SRCEXTERN)/gsl*.$(ARCHIVE) | grep -Eo '[0-9]+(.[0-9]+)*')
DEP_GSL_BUILD := $(abspath $(BUILDEXTERN)/gsl-$(DEP_GSL_V))
DEP_GSL_SRC := $(SRCEXTERN)/gsl-$(DEP_GSL_V)

# Boost
DEP_BOOST_V=$(shell echo $(SRCEXTERN)/boost*.$(ARCHIVE_BOOST) | grep -Eo '[0-9]+(.[0-9]+)+')
DEP_BOOST_BUILD := $(abspath $(BUILDEXTERN)/boost_$(DEP_BOOST_V))
DEP_BOOST_SRC := $(SRCEXTERN)/boost_$(DEP_BOOST_V)

BOOST_NOBUILD=atomic,chrono,container,context,coroutine,date_time,exception,filesystem,graph,graph_parallel,iostreams,locale,log,math,mpi,program_options,python,random,regex,serialization,signals,system,test,thread,timer,wave

CC=g++
CFLAGS=-c -Wall
SPECTRAST=spectrast


LDFLAGS= -lm \
		 -lpthread \
		 -L$(DEP_EXPAT_BUILD)/lib -lexpat \
		 -L$(DEP_GSL_BUILD)/lib  -lgsl -lgslcblas\
		 -L$(DEP_ZLIB_BUILD)/lib -lz

IFLAGS= -I/usr/include \
	   	-I$(DEP_EXPAT_BUILD)/include  \
		-I$(DEP_ZLIB_BUILD)/include \
	   	-I$(DEP_GSL_BUILD)/include \
	   	-I$(DEP_BOOST_BUILD)/include

# SpectraST sources
CPP_FILES_SPECTRAST := $(wildcard src/spectrast/*.cpp)
OBJ_FILES_SPECTRAST := $(addprefix obj/spectrast/,$(notdir $(CPP_FILES_SPECTRAST:.cpp=.o))) 

# Ramp sources
CPP_FILES_RAMP := $(wildcard src/ramp/*.cpp)
OBJ_FILES_RAMP := $(addprefix obj/ramp/,$(notdir $(CPP_FILES_RAMP:.cpp=.o)))


all: $(DEP_EXPAT_BUILD) $(DEP_ZLIB_BUILD) $(DEP_GSL_BUILD) $(DEP_BOOST_BUILD) $(SPECTRAST)

$(DEP_EXPAT_BUILD): 
	mkdir -p $(DEP_EXPAT_BUILD)
	cd $(SRCEXTERN) ; tar xf expat*$(ARCHIVE) 
	cd $(DEP_EXPAT_SRC);  ./configure --prefix=$(DEP_EXPAT_BUILD); make; make install

$(DEP_ZLIB_BUILD):
	mkdir -p $(DEP_ZLIB_BUILD)
	cd $(SRCEXTERN) ; tar xf zlib*$(ARCHIVE)
	cd $(DEP_ZLIB_SRC);  ./configure --prefix=$(DEP_ZLIB_BUILD); make; make install

$(DEP_GSL_BUILD):
	mkdir -p $(DEP_GSL_BUILD)
	cd $(SRCEXTERN) ; tar xf gsl*$(ARCHIVE)
	cd $(DEP_GSL_SRC);  ./configure --prefix=$(DEP_GSL_BUILD) --with-pic --disable-shared; make; make install

$(DEP_BOOST_BUILD):
	mkdir -p $(DEP_BOOST_BUILD)
	cd $(SRCEXTERN) ; tar xf boost*$(ARCHIVE_BOOST)
	cd $(DEP_BOOST_SRC);  ./bootstrap.sh --prefix=$(DEP_BOOST_BUILD) --without-libraries=$(BOOST_NOBUILD); ./b2 install

$(SPECTRAST): $(OBJ_FILES_RAMP) $(OBJ_FILES_SPECTRAST) 
	    $(CC) $(IFLAGS) -o $@ $^ $(LDFLAGS) -DSTANDALONE_LINUX -D_FILE_OFFSET_BITS=64 

obj/spectrast/%.o: src/spectrast/%.cpp
	mkdir -p obj/spectrast
	$(CC) $(IFLAGS) $(CFLAGS) -o $@ $<  -DSTANDALONE_LINUX -D_FILE_OFFSET_BITS=64


obj/ramp/%.o: src/ramp/%.cpp
	mkdir -p obj/ramp
	$(CC) $(IFLAGS) $(CFLAGS) -o $@ $<  -DSTANDALONE_LINUX -D_FILE_OFFSET_BITS=64


.PHONY: clean
clean:
	rm -f $(SPECTRAST)
	rm -rf $(OBJ)
	rm -rf $(BUILD)
	rm -rf $(DEP_EXPAT_SRC)
	rm -rf $(DEP_ZLIB_SRC)
	rm -rf $(DEP_GSL_SRC)
	rm -rf $(DEP_BOOST_SRC)
