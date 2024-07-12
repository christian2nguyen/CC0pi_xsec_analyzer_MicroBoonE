CXX = g++
CXXFLAGS = -g -Wall -fPIC -Wno-unused-variable
ROOTFLAGS = `root-config --cflags --glibs --libs` -lTreePlayer -lEG -lMinuit
XGBOOST=/exp/uboone/app/users/mastbaum/xgboost/build
CC0PI_ANALYZER_FLAGS = -I${XGBOOST}/include -L${XGBOOST}/lib64 -lxgboost -g


# make a binary for every .cxx file
 all : $(patsubst %.cpp, %.o, $(wildcard *.cpp))  test_unfolding_executable FlatTree_2DAnaylizer GeneratorOverlay_executable
# cc0pi_analyzer_org
# # rule for each one
# 
%: %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $(CC0PI_ANALYZER_FLAGS) $@ $< ../includes/*.o
	
 %.o : %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $*.o  -c $*.cpp 
	$(CXX( $(CXXFLAGS) $(ROOTFLAGS) -o $* $*.o ../includes/*.o 

stv_root_dict.o:
	$(RM) stv_root_dict*.*
	rootcling -f stv_root_dict.cc -c LinkDef.h
	$(CXX) $(shell root-config --cflags --libs) -O3 \
	-fPIC -o stv_root_dict.o -c stv_root_dict.cc
	$(RM) stv_root_dict.cc
	
test_unfolding_executable: test_unfolding_executable.cpp
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(CC0PI_ANALYZER_FLAGS) -O3 -o $@ $^ ../includes/*.o
	
FlatTree_2DAnaylizer: FlatTree_2DAnaylizer.cpp	
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(CC0PI_ANALYZER_FLAGS) -O3 -o $@ $^ ../includes/*.o
	
GeneratorOverlay_executable: GeneratorOverlay_executable.cpp 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(CC0PI_ANALYZER_FLAGS) -O3 -o $@ $^ ../includes/*.o
	
	
.PHONY: clean

.INTERMEDIATE: stv_root_dict.o

clean:
	rm -f $(wildcard *.o) $(patsubst %.cpp, %, $(wildcard *.cpp)) test_unfolding_executable FlatTree_2DAnaylizer cc0pi_analyzer GeneratorOverlay_executable 
