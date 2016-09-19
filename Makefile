# make file to compile and link a test program that uses the verbmenu library

CCFILE       = shelfmc_stripped.cc  functions.cc 

# Define filename suffixes
ObjSuf        = .o
SrcSuf        = .cc
IncSuf        = .hh
ExeSuf        = .exe
DllSuf        = .so
OutPutOpt     = -o

# Define the compile and link commands
CXX           = g++ 
CXXFLAGS      = -c -g -O2 -o $(OBJ) -Wall -fPIC -I`${ROOTSYS}/bin/root-config --incdir`
LD            = g++ 
LDFLAGS       = -O 

# Define root libraries
#ROOTLIBS      = -L$(ROOTSYS)/lib 
ROOTLIBS      = `${ROOTSYS}/bin/root-config --libs` 

# Define all my libraries	
LIBS          = $(ROOTLIBS) 

# Define shortcuts for compiling and linking
COMPILE = $(CXX) $(CXXFLAGS) $(SCRATCH)
LINK = 	$(LD) $(LDFLAGS) $(OBJ) $(LIBS) $(OutPutOpt) $(EXE)
#------------------------------------------------------------------------------

# Define the file names
OBJ       = shelfmc_stripped.obj
SRC       = $(CCFILE)
EXE       = shelfmc_stripped.exe
SCRATCH   = tmp_stripped.cc
INC       = shelfmc_stripped.inc 

all:            $(EXE)

# Update the executable if the object file has changed
$(EXE): 	$(OBJ)
		$(LINK)


# Update the object file if the source, or include file changed
$(OBJ):     $(SRC)
	      cat $(CCFILE) > $(SCRATCH)
		$(COMPILE)

clean:
	@rm -f $(OBJ) $(EXE) core* shelfmc_stripped.exe shelfmc_stripped.obj shelfmc_stripped.root $(SCRATCH)

.SUFFIXES: .$(SrcSuf)

.$(SrcSuf).$(ObjSuf):	
	$(CXX) $(CXXFLAGS) -c $<	







