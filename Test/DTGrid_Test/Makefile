LIBS = 
CXXFLAGS += -O3 -I../..
LDFLAGS =

SRC = SourceFiles/DTGrid_Test.cpp
OBJ = $(SRC:%.cpp=%.o)

all:	$(OBJ)
	$(CXX) $(LDFLAGS) -o DTGrid_Test.exe $(OBJ) $(LIBS) 

install:	all
	cp DTGrid_Test.exe $(DTGridBinariesPath)

clean:
	rm -f SourceFiles/*.o DTGrid_Test.exe

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@


