# -- TARGET --------------------------------------------------------------------
#
# Determine details about what platform we are targeting the build for. The 
# following variables are filled in based on the target:
#
# ARCH      The build machine architecture targeted (i686, x86_64)
# VENDOR    The distributor of the kernel software and package management 
#           system the build is targeted for (e.g. redhat, ubuntu, mingw32, w64)
# SYSTEM    The target kernel/system (linux, darwin, mingw32)
# OS        The host operating system TPP is being built on. One of Linux, 
#           Windows_NT, or Darwin
#
# There are a number of ways for deducing this depending on the system you are
# on and what you are building for. For TPP's usual platforms here are the 
# common ways and values reported:
#
#                    Red Hat Enterprise  MinGW/MinGW-64      OS X/Darwin
# ${OS}                                  Windows_NT          ?
# uname              Linux               MINGW32_NT-6.1      Darwin
# uname -s           Linux               MINGW32_NT-6.1      Darwin
# uname -o           GNU/Linux           Msys                Darwin
# uname -m           x86_64              i686                x86_64
# uname -p           x86_64              unknown             i386
# uname -i           x86_64              unknown             <error>
# gcc -dumpmachine   x86_64-redhat-linux mingw32 |           i686-apple-darwin11
#                                        i686-w64-mingw32 |
#                                        x86_64-w64-mingw32
#
# We'll use the current platform's <arch>-<vendor>-<os> triplet as reported by
# the gcc compiler for the build information.  (And of course msys/mingw32 
# doesn't report it as a triplet.)
#
# C++ STANDARD 

CXXSTD := 
GCCVERSION := $(shell $(CXX) -dumpversion)

ifeq ("$(GCCVERSION)","5.2.0")
CXXSTD := -std=c++11
endif

TARGET := $(shell $(CXX) -dumpmachine)
ifeq ($(TARGET),mingw32)
   $(error Error: Should not build using MinGW gcc tool chain instead use MingGW-w64)
endif

# Fill in the details
ifneq (,$(findstring linux, $(TARGET)))
   OS := Linux
else ifneq (,$(findstring darwin, $(TARGET)))
   OS := Darwin
else ifneq (,$(findstring mingw, $(TARGET)))
   OS := $(OS)
else ifneq (,$(findstring msys, $(TARGET)))
   OS := $(OS)
else
   $(error Unable to determine target platform using $$(CXX) -dumpmachine)
endif

# -- COMPILER DETAILS  ---------------------------------------------------------
#
# Compiler and common compiler flags

CXX := g++ -static -static-libgcc -static-libstdc++
CXXFLAGS ?=

#GDB = 1 #uncomment to compile for gdb debugger
#DEBUG = 1 #uncomment to set DEBUG flag
# Debug or optimize?

CXXFLAGS += $(if $(GDB),-g,-O2)
CXXFLAGS += $(if $(DEBUG),-DDEBUG)

CXXFLAGS += $(CXXSTD)

# Turn all warnings into errors but ignore deprecated function calls for now
# as they are treated as errors and can't be be ignored by the 
# option -Wno-errors=deprecated.
#
# TODO: remove the deprecated calls
CXXFLAGS += -Wno-deprecated

# Always include large file support
CXXFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE

# Make a LGPL instead of GPL version of TPP (deprecated)
CXXFLAGS += $(if $(LGPL_SUPPORT),-D__LGPL__)

# Include MZ5 support in mzParser (requires HDF5 libraries)
CXXFLAGS += $(if $(MZ5_SUPPORT),-DMZP_MZ5)

# TODO: get rid of this if possible
# CXXFLAGS += -DTPPLIB

CXXFLAGS += -DSTANDALONE_LINUX
# System specific
ifeq ($(OS),Linux)
   CXXFLAGS += -D__LINUX__
   # compile with position independent code to allow for library inclusion
   # CXXFLAGS += -fPIC
endif
ifeq ($(OS),MingW)
   CXXFLAGS += -D__MINGW__ -D_USE_32BIT_TIME_T -D_GLIBCXX_USE_WCHAR_T
endif
ifeq ($(OS),Windows_NT)
   CXXFLAGS += -DWIN32 -D_USE_MATH_DEFINES 
endif
ifeq ($(OS),Darwin)
   CXXFLAGS += -D__LINUX__ -ftemplate-depth=256
endif


# -- LINKER DETAILS  -----------------------------------------------------------
#
# Linker and common linker flags

LD = g++ -static -static-libgcc -static-libstdc++
LDFLAGS += $(if $(DEBUG)||$(GDB),,-s)     # Strip executables
LDFLAGS += -static -static-libgcc -static-libstdc++


