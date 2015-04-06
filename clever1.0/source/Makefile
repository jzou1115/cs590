# CE method
# written by Seungsoo Hahn

CC = g++
CXXFLAGS = -O2 -static -D__LINUX -Wno-deprecated
CPPFLAGS =
DEPEND = makedepend $(CXXFLAGS)

SOURCES =   argparser.cpp    ccfs.cpp cemethod.cpp      ctset.cpp\
               ccfdef.cpp   cctrl.cpp   cinput.cpp coperators.cpp\
            exception.cpp ccflist.cpp  mcsimul.cpp    cstring.cpp\
                cutil.cpp

OBJS := $(SOURCES:.cpp=.o) 


all: CEEnergy CETrFile GenSeqs

depend: $(SOURCES)
	$(DEPEND) $(SOURCES)

OBJSE := Main_CEEnergy.o

OBJSF := Main_CETrFile.o

OBJSG := Main_GenSeqs.o

CEEnergy: $(OBJS) $(OBJSE) 
	$(CC) $(CXXFLAGS) -o $@ $(OBJS) $(OBJSE)

CETrFile: $(OBJS) $(OBJSF)
	$(CC) $(CXXFLAGS) -o $@ $(OBJS) $(OBJSF)

GenSeqs: $(OBJS) $(OBJSG)
	$(CC) $(CXXFLAGS) -o $@ $(OBJS) $(OBJSG)

clean:
	- $(RM) CEEnergy CETrFile GenSeqs *.o

