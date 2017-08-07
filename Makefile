INCDIR=include
SRCDIR=src
OBJDIR=obj
EXE=sean

CC=g++
CFLAGS=-g -Wall -Wconversion -I$(INCDIR)
ROOTFLAGS=`root-config --cflags --glibs` -lMathMore

_DEPS = Experiment.h Config.h CrossSection.h Target.h Absorption.h
DEPS = $(patsubst %,$(INCDIR)/%,$(_DEPS))

_OBJ = Experiment.o main.o CrossSection.o Target.o Absorption.o
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))
#
#$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
#	$(CC) -c -o $@ $< $(CFLAGS) $(ROOTFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(ROOTFLAGS) 

$(EXE): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(ROOTFLAGS)

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(EXE)
