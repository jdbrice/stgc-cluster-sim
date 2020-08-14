
SIM_OBJS = stgc.o
SIM_EXE = sim

FIN_OBJS = finder.o
FIN_EXE = finder

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMathMore
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++
FLAGS = -std=c++11 -fno-inline -Wall -g $(INCFLAGS) $(LDFLAGS) 

COMPILE = $(CXX) $(FLAGS) -c 

all: $(SIM_EXE) $(FIN_EXE)

$(SIM_EXE): $(SIM_OBJS)
	$(CXX) -fno-inline -Wl -lpthread -o $(SIM_EXE) $(SIM_OBJS) $(ROOTFLAGS) $(ROOTLIBS)
$(FIN_EXE): $(FIN_OBJS)
	$(CXX) -fno-inline -Wl -lpthread -o $(FIN_EXE) $(FIN_OBJS) $(ROOTFLAGS) $(ROOTLIBS)


%.o: %.cpp
	$(COMPILE)  $< 

