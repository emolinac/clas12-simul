ifndef ROOTSYS
$(error "Please set the variable ROOTSYS")
endif

ifndef CERN_LIB
$(error "Please set the variable CERN_LIB")
endif

ifndef CLAS_PACK
$(error "Please set the variable CLAS_PACK")
endif

ARCH = $(shell root-config --ldflags)

MKDIR_P := mkdir -p

OBJDIR := ./obj
BINDIR := ./bin
SRCFDIR := ./src_f
SRCCDIR := ./src

EXE := $(BINDIR)/lepto.exe

CXX := g++
FF := gfortran

FFLAGS =  $(ARCH) -g -O -fno-automatic -fno-second-underscore
FFLAGS += -ffixed-line-length-none -funroll-loops -Wunused
FFLAGS += -I$(CLAS_PACK)/include

CCFLAGS =  $(ARCH) -g -Wall -fPIC -Wno-deprecated
CCFLAGS += $(shell root-config --cflags)
CCFLAGS += -I$(ROOTSYS)/include -I$(CLAS_PACK)/include

ROOTLIB := $(shell root-config --glibs)
CERNLIB	:= -L$(CERN_LIB) -Wl,-static -lmathlib -lpacklib -lkernlib -Wl,-dy -lm -lnsl -lcrypt -ldl -lgfortran
CLASLIB := -L$(CLAS_LIB) -lbosio
MYLIBS := -Wl,-rpath, -lm -lgfortran
LIBS := $(ROOTLIB) $(CERNLIB) $(CLASLIB) $(MYLIBS)

FSRCS = $(wildcard $(SRCFDIR)/*.f)
FOBJS := $(FSRCS:%.f=%.o)
FOBJS := $(FOBJS:$(SRCFDIR)/%=$(OBJDIR)/%)

CCSRCS = $(wildcard $(SRCCDIR)/*.cc)
CCOBJS := $(CCSRCS:%.cc=%.o)
CCOBJS := $(CCOBJS:$(SRCCDIR)/%=$(OBJDIR)/%)

all: $(EXE)

$(EXE): $(CCOBJS) $(FOBJS)
	$(MKDIR_P) $(BINDIR)
	$(CXX) $(CCFLAGS) -o $@ $^ $(LIBS)

$(OBJDIR)/%.o: $(SRCCDIR)/%.cc
	$(MKDIR_P) $(OBJDIR)
	@echo "Compiling " $@
	$(CXX) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCFDIR)/%.f
	$(MKDIR_P) $(OBJDIR)
	@echo "Compiling " $@
	$(FF) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(EXE) $(CCOBJS) $(FOBJS)
