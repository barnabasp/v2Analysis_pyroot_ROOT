# At zeroth order, the user does not have to below this line.
WrkDir := $(shell pwd)
ExeDir = exe
LibDir = lib
ObjDir = $(ExeDir)/object
DepDir = $(ObjDir)/dependencies

CXX     = g++
LD      = g++
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit
CFLAGS  = -O0 -g -Wall -fPIC -fno-inline $(ROOTCFLAGS) -I$(WrkDir)/$(RdrDir)
LDFLAGS = -O0 # beware -m32 may have to be deleted/added
SYSLIBS = 

COMMON_SOURCES = particle_tree.C
COMMON_OBJECTS = $(addprefix $(ObjDir)/,$(addsuffix .o,$(notdir $(basename $(COMMON_SOURCES)))))
SOURCES = $(COMMON_SOURCES)
ALL_SOURCES = $(sort $(SOURCES))

# Check system for shared library extension
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    SHARED_LIB_EXT := dylib
else
    SHARED_LIB_EXT := so
endif

define COMPILE_TEMPLATE
-include $(DepDir)/$(notdir $(basename $(1))).d
$(ObjDir)/$(notdir $(basename $(1))).o:
	@echo ""
	$(CXX) $(CFLAGS) -c -MD -MP -MF $(DepDir)/$(notdir $(basename $(1))).d $(1) -o $$@
endef
$(foreach source, $(ALL_SOURCES), $(eval $(call COMPILE_TEMPLATE,$(source))))

$(LibDir)/libcommon.$(SHARED_LIB_EXT): $(COMMON_OBJECTS)
	@echo ""
	$(LD) -shared -o $@ $^ $(ROOTLIBS) $(SYSLIBS)

all: $(LibDir)/libcommon.$(SHARED_LIB_EXT)

clean:
	@rm -f $(LibDir)/*.$(SHARED_LIB_EXT) $(ObjDir)/*.o $(DepDir)/*.d

