MKDIR_P := mkdir -p

OBJDIR := ./obj
BINDIR := ./bin
SRCCDIR := ./src

EXE := $(BINDIR)/txt2part.exe

CC := gcc

CCFLAGS := -m64 -Wall -fPIC -Wno-deprecated -g -O0

CERN_LIBRARIES := -L$(CERN_LIB) -Wl,-static -lmathlib -lpacklib -lkernlib -Wl,-dy -lm -lnsl -lcrypt -ldl -lgfortran
CLAS_LIBRARIES := -L$(CLAS_LIB) -lscaler -lclasutil -lc_bos_io -lmapmanager -lfputil -lbos -lfpack -lc_cern
CCLIBS := $(CLAS_LIBRARIES) $(CERN_LIBRARIES)

CLAS_INCLUDES := -I$(CLAS_PACK)/include
INCLUDES := $(CLAS_INCLUDES)

CCSRCS = $(wildcard $(SRCCDIR)/*.c)
CCOBJS := $(CCSRCS:%.c=%.o)
CCOBJS := $(CCOBJS:$(SRCCDIR)/%=$(OBJDIR)/%)

all: $(EXE)

$(EXE): $(CCOBJS)
	$(MKDIR_P) $(BINDIR)
	@echo "Linking" $@
	$(CC) $(CCFLAGS) $(INCLUDES) -o $@ $^ $(CCLIBS)

$(OBJDIR)/%.o: $(SRCCDIR)/%.c
	$(MKDIR_P) $(OBJDIR)
	@echo "Compiling" $@
	$(CC) $(CCFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(EXE) $(CCOBJS)
