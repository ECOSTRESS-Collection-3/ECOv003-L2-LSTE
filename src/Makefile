# !gmake
# DEBUG = 1
#
# !Makefile Name: Makefile
#
# !Description: Makefile for L2_PGE Land Surface Temperature and Emissivity.
#               This makefile used to generate l2_tes_pge
#
# !Env Variables:
#
#      ENV VARIABLES DESCRIPTION
#      ~~~~~~~~~~~~~ ~~~~~~~~~~~
#      CC            The c compiler
#      CFLAGS        Compiler flags
#      HDFINC        Include directory of HDF
#      HDF5INC       Include directory of HDF5
#      HDFLIB        Library for HDF
#      HDF5LIB       Library for HDF5
#
# !Revision History:
#  Date:        PR#     Author          Description
#  -------      -----   ------          -----------
#  05/03/16             R. Freepartner  Initial Version
###################################################################
# Auto-detect conda env when CONDA_PREFIX is not exported (common in non-activated shells)
DEFAULT_ENV_NAME ?= ECOv003-L2-LSTE
CONDA_BASE_DETECTED := $(shell command -v conda >/dev/null 2>&1 && conda info --base 2>/dev/null || echo "")
CONDA_PREFIX_DETECTED := $(shell \
  if [ -n "$(CONDA_BASE_DETECTED)" ] && [ -d "$(CONDA_BASE_DETECTED)/envs/$(DEFAULT_ENV_NAME)" ]; then echo "$(CONDA_BASE_DETECTED)/envs/$(DEFAULT_ENV_NAME)"; \
  elif [ -n "$(CONDA_PREFIX)" ]; then echo "$(CONDA_PREFIX)"; \
  elif [ -d "/opt/homebrew/Caskroom/miniforge/base/envs/$(DEFAULT_ENV_NAME)" ]; then echo "/opt/homebrew/Caskroom/miniforge/base/envs/$(DEFAULT_ENV_NAME)"; \
  fi)

ifneq ($(CONDA_PREFIX_DETECTED),)
CONDA_PREFIX := $(CONDA_PREFIX_DETECTED)
export CONDA_PREFIX
PKG_CONFIG_PATH := $(CONDA_PREFIX)/lib/pkgconfig:$(PKG_CONFIG_PATH)
export PKG_CONFIG_PATH
endif

# Auto-discover dependencies using pkg-config
# This makes the build portable across any system with installed packages

PKG_CHECK = $(shell pkg-config --cflags $(1) 2>/dev/null)
PKG_LIBS = $(shell pkg-config --libs $(1) 2>/dev/null)

# Skip all dependency checks when running targets that don't need them
CLEANING := $(filter clean clean_o clean_exe install,$(MAKECMDGOALS))
ifneq ($(CLEANING),)
HDFINC  :=
HDFLIB  :=
HDF5INC :=
HDF5LIB :=
XMLINC  :=
XMLLIB  :=
EXPATLIB :=
GRIBINC :=
GRIBLIB :=
else

# HDF4 (required) - provides hdf.h, mfhdf.h, hdfi.h
# Try pkg-config first, then check conda/homebrew paths (HDF4 often doesn't have .pc file)
HDF4_PREFIX := $(shell if [ -n "$$CONDA_PREFIX" ]; then echo $$CONDA_PREFIX; fi)
ifeq ($(HDF4_PREFIX),)
  HDF4_PREFIX := $(shell command -v brew >/dev/null 2>&1 && brew --prefix hdf4 2>/dev/null || echo "")
endif

ifeq ($(HDF4_PREFIX),)
  HDFINC := $(call PKG_CHECK,hdf)
  HDFLIB := $(call PKG_LIBS,hdf)
else
  HDFINC := -I$(HDF4_PREFIX)/include
  HDFLIB := -L$(HDF4_PREFIX)/lib -lmfhdf -ldf
endif

ifeq ($(HDFINC),)
$(error HDF4 not found. Use: conda create -n ecostress -c conda-forge hdf4 && conda activate ecostress)
endif

# HDF5 (required) - provides hdf5.h
HDF5INC := $(call PKG_CHECK,hdf5)
HDF5LIB := $(call PKG_LIBS,hdf5)
ifeq ($(HDF5INC),)
$(error HDF5 not found. Install with: brew install hdf5 (macOS) or apt-get install libhdf5-dev (Ubuntu))
endif

# LibXML2 (required) - provides libxml/xmlwriter.h
# LibXML2 often doesn't provide .pc file in conda, so check conda/homebrew paths
XMLINC := $(call PKG_CHECK,libxml-2.0)
ifeq ($(XMLINC),)
  # Try conda/homebrew paths
  ifdef CONDA_PREFIX
    XMLINC := -I$(CONDA_PREFIX)/include/libxml2
    XMLLIB := -L$(CONDA_PREFIX)/lib -lxml2
  else
    XMLINC := $(shell if command -v brew >/dev/null 2>&1; then echo -I$$(brew --prefix libxml2)/include/libxml2; fi)
    XMLLIB := $(shell if command -v brew >/dev/null 2>&1; then echo -L$$(brew --prefix libxml2)/lib -lxml2; fi)
  endif
endif
ifeq ($(XMLINC),)
$(error libxml2 not found. Use: conda create -n ecostress -c conda-forge hdf4 hdf5 libxml2 eccodes pkg-config && conda activate ecostress)
endif
XMLLIB_VAL := $(XMLLIB)
ifeq ($(XMLLIB_VAL),)
  XMLLIB_VAL := $(call PKG_LIBS,libxml2-2.0)
endif
ifeq ($(XMLLIB_VAL),)
  XMLLIB_VAL := -lxml2
endif
XMLLIB := $(XMLLIB_VAL)

# Expat (required) - provides expat.h
# Expat often doesn't provide .pc file in conda, so check conda/homebrew paths
EXPATLIB := $(call PKG_LIBS,expat)
ifeq ($(EXPATLIB),)
  # Try conda/homebrew paths
  ifdef CONDA_PREFIX
    EXPATLIB := -L$(CONDA_PREFIX)/lib -lexpat
  else
    EXPATLIB := $(shell if command -v brew >/dev/null 2>&1; then echo -L$$(brew --prefix expat)/lib -lexpat; fi)
  endif
endif
ifeq ($(EXPATLIB),)
$(error expat not found. Use: conda create -n ecostress -c conda-forge hdf4 hdf5 libxml2 eccodes pkg-config && conda activate ecostress)
endif

# GRIB API / ecCodes (required) - provides grib_api.h
# Try ecCodes first (modern), fall back to grib_api (legacy)
GRIBINC := $(call PKG_CHECK,eccodes)
GRIBLIB := $(call PKG_LIBS,eccodes)
ifeq ($(GRIBLIB),)
  GRIBINC := $(call PKG_CHECK,grib_api)
  GRIBLIB := $(call PKG_LIBS,grib_api)
endif
ifeq ($(GRIBLIB),)
  # Fallback for conda setups where pkg-config lookup is unreliable
  ifdef CONDA_PREFIX
    ifneq ($(wildcard $(CONDA_PREFIX)/include/eccodes.h),)
      GRIBINC := -I$(CONDA_PREFIX)/include
      GRIBLIB := -L$(CONDA_PREFIX)/lib -leccodes
    endif
  endif
endif
ifeq ($(GRIBLIB),)
$(error ecCodes/GRIB API not found. Install with: brew install eccodes (macOS) or apt-get install libeccodes-dev (Ubuntu))
endif

endif # ifneq ($(CLEANING),)

##################
# EXPAT settings
# Adds support for checking lexical aspects of namespaces
XP_NS = -DXML_NS
# Adds support for processing DTDs
XP_DTD = -DXML_DTD
# Byte order macro
#This is for little endian machines
XP_BO = -DXML_BYTE_ORDER=12

XP_DEFINES = $(XP_NS) $(XP_DTD) $(XP_BO) $(XP_UNI) $(XP_UNI_WC) $(XP_MM)
##################

ifdef DEBUG
ADD_CFLAGS = -g -Wall
else
ADD_CFLAGS = -O3 -Wall
endif

BUILDDIR = ../build
BINDIR  = ../bin

# Object files (built into ../build/)
TES_OBJ = \
    $(BUILDDIR)/tes_main.o \
    $(BUILDDIR)/asterged.o \
    $(BUILDDIR)/cloud.o \
    $(BUILDDIR)/config.o \
    $(BUILDDIR)/describe.o \
    $(BUILDDIR)/error.o \
    $(BUILDDIR)/fileio.o \
    $(BUILDDIR)/interps.o \
    $(BUILDDIR)/lste_lib.o \
    $(BUILDDIR)/matrix.o \
    $(BUILDDIR)/maptogrid.o \
    $(BUILDDIR)/metadata.o \
    $(BUILDDIR)/rttov_util.o \
    $(BUILDDIR)/smooth2d.o \
    $(BUILDDIR)/tes_util.o \
    $(BUILDDIR)/tg_wvs.o

ALLOBJ = $(TES_OBJ)

# Combine the includes using pkg-config results
INC = -I. -I../include $(HDFINC) $(HDF5INC) $(XMLINC) $(GRIBINC) -I/usr/include

# Header files (now in ../include/)
HDR = ../include/asterged.h \
      ../include/cloud.h \
      ../include/config.h \
      ../include/describe.h \
      ../include/fileio.h \
      ../include/interps.h \
      ../include/lste_lib.h \
      ../include/maptogrid.h \
      ../include/matrix.h \
      ../include/metadata.h \
      ../include/rttov_util.h \
      ../include/smooth2d.h \
      ../include/tes_util.h \
      ../include/tg_wvs.h
      
# Combine the libraries using pkg-config results
LIB = $(HDFLIB) $(HDF5LIB) $(GRIBLIB) $(XMLLIB) $(EXPATLIB) -ldl -ldf -ljpeg -lz -lm 

# Executable goes to ../bin/
TARGET = $(BINDIR)/L2_PGE

all: $(BUILDDIR) $(BINDIR) $(TARGET)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(TARGET): $(TES_OBJ) $(HDR)
	$(CC) $(ADD_CFLAGS) $(XP_DEFINES) $(TES_OBJ) $(LIB) -o $(TARGET)

$(BUILDDIR)/%.o: %.c $(HDR)
	$(CC) $(CFLAGS) $(ADD_CFLAGS) $(INC) -c $< -o $@

# Delete object files
clean: clean_exe clean_o

clean_o:
	rm -f $(ALLOBJ)

# Delete the executable
clean_exe:
	rm -f $(TARGET)

# Install the executable (default prefix: /usr/local, override with PREFIX=...)
PREFIX ?= /usr/local
install:
	@if [ ! -f "$(TARGET)" ]; then \
		echo "Error: $(TARGET) not found. Run 'make' first."; \
		exit 1; \
	fi
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $(TARGET) $(DESTDIR)$(PREFIX)/bin/

###################################################################
