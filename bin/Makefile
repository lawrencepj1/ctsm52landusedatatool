all: ctsm52landusedatatool

# Check for the netcdf library and include directories 
ifeq ($(LIB_NETCDF),$(null))
  LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
  INC_NETCDF := /usr/local/include
endif

ifeq ($(MOD_NETCDF),$(null))
  MOD_NETCDF := $(LIB_NETCDF)
endif

ctsm52landusedatatool: ../src/ctsm52landusedatatool.c
	icc -o ctsm52landusedatatool ../src/ctsm52landusedatatool.c -lm -mcmodel=medium -lnetcdf

#	cc -o ctsm52landusedatatool ../src/ctsm52landusedatatool.c -lm -mcmodel=medium -L$(LIB_NETCDF) -I$(INC_NETCDF) -lnetcdf
