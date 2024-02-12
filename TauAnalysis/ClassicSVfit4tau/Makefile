CXX      = g++
FC       = gfortran

WFLAGS   = -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -Wno-vla
CXXFLAGS = $(WFLAGS) -fvisibility-inlines-hidden -fno-math-errno --param vect-max-version-for-alias-checks=50 \
           -Xassembler --compress-debug-sections -fno-crossjumping -msse3 -felide-constructors -fmessage-length=0 \
           -fdiagnostics-show-option -pipe -fPIC $(shell root-config --cflags)

FCWFLAGS = -Wunused -Wuninitialized
FCFLAGS  = $(FCWFLAGS) -c -fno-second-underscore -cpp -fPIC

LDFLAGS  = -flto

ifdef DEBUG
  CXXFLAGS += -Og -ggdb3
  FCFLAGS  += -Og
else
  CXXFLAGS += -O3 -ggdb0
  LDFLAGS  += -s
  FCFLAGS  += -O3
endif

BASEDIR = TauAnalysis/ClassicSVfit4tau
VAMPDIR = VAMP/install/include/vamp

SOURCE_PATH = $(BASEDIR)/src
BIN_PATH    = $(BASEDIR)/bin

OBJ_PATH  = $(BASEDIR)/obj
DEP_PATH  = $(BASEDIR)/dep
LIB_PATH  = $(BASEDIR)/lib
EXEC_PATH = $(BASEDIR)/exec

SRC_EXT = cc
OBJ_EXT = o
DEP_EXT = d
LIB_EXT = so
FC_EXT  = f

INCLUDES = -I. -I$(VAMPDIR)
LIBS     = $(shell root-config --libs) -lGenVector -lTauAnalysis_ClassicSVfit -lvamp

SRCS          = $(shell find $(SOURCE_PATH) -name '*.$(SRC_EXT)' -exec basename {} \;)
OBJS          = $(SRCS:%.$(SRC_EXT)=$(OBJ_PATH)/%.$(OBJ_EXT))
DEPS          = $(SRCS:%.$(SRC_EXT)=$(DEP_PATH)/%.$(DEP_EXT))
TRGT_SRCS     = $(shell find $(BIN_PATH) -name '*.$(SRC_EXT)' -exec basename {} \;)
TRGT_OBJS     = $(TRGT_SRCS:%.$(SRC_EXT)=$(OBJ_PATH)/%.$(OBJ_EXT))
TRGT          = $(TRGT_SRCS:%.$(SRC_EXT)=$(EXEC_PATH)/%)
TRGT_LIB_BASE = $(subst /,_,$(BASEDIR))
TRGT_LIB_PATH = $(LIB_PATH)/lib$(TRGT_LIB_BASE).$(LIB_EXT)
SRCS_FC       = $(shell find $(SOURCE_PATH) -name '*.$(FC_EXT)' -exec basename {} \;)
OBJS_FC       = $(SRCS_FC:%.$(FC_EXT)=$(OBJ_PATH)/%.$(OBJ_EXT))

all: $(TRGT)

$(TRGT): $(EXEC_PATH)/%: $(OBJ_PATH)/%.$(OBJ_EXT) $(TRGT_LIB_PATH)
	@mkdir -p $(EXEC_PATH)
	@$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS) -l$(TRGT_LIB_BASE)

$(TRGT_LIB_PATH): $(OBJS) $(OBJS_FC)
	@mkdir -p $(LIB_PATH)
	@$(CXX) $(CXXFLAGS) -shared $^ -o $@ $(LDFLAGS) $(LIBS)

$(OBJ_PATH)/%.$(OBJ_EXT): $(SOURCE_PATH)/%.$(SRC_EXT)
	@mkdir -p $(@D)
	@mkdir -p $(DEP_PATH)
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MF $(patsubst $(OBJ_PATH)/%.$(OBJ_EXT),$(DEP_PATH)/%.$(DEP_EXT),$@) -c $< -o $@

$(OBJ_PATH)/%.$(OBJ_EXT): $(BIN_PATH)/%.$(SRC_EXT)
	@mkdir -p $(@D)
	@mkdir -p $(DEP_PATH)
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MF $(patsubst $(OBJ_PATH)/%.$(OBJ_EXT),$(DEP_PATH)/%.$(DEP_EXT),$@) -c $< -o $@

$(OBJ_PATH)/%.$(OBJ_EXT): $(SOURCE_PATH)/%.$(FC_EXT)
	@mkdir -p $(@D)
	@mkdir -p $(DEP_PATH)
	@$(FC) $(FCFLAGS) $(INCLUDES) -MMD -MF $(patsubst $(OBJ_PATH)/%.$(OBJ_EXT),$(DEP_PATH)/%.$(DEP_EXT),$@) -c $< -o $@

.PHONY: clean

clean:
	@rm -rf $(OBJ_PATH) $(DEP_PATH) $(EXEC_PATH)
