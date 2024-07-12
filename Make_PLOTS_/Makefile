CXX = g++
CXXFLAGS = -g -Wall -fPIC -Wno-unused-variable
ROOTFLAGS = `root-config --cflags --glibs --libs` -lTreePlayer -lEG -lMinuit
XGBOOST=/exp/uboone/app/users/mastbaum/xgboost/build
CC0PI_ANALYZER_FLAGS = -I${XGBOOST}/include -L${XGBOOST}/lib64 -lxgboost -g


# make a binary for every .cxx file
 all : $(patsubst %.cpp, %.o, $(wildcard *.cpp)) Make_Plots Make_Plots_2D tutorial_slice_plots_test Make_Plots_2D_slicing_test
# cc0pi_analyzer_org
# # rule for each one
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
	
Make_Plots_2D_slicing_test:Make_Plots_2D_slicing_test.cpp 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(CC0PI_ANALYZER_FLAGS) -O3 -o $@ $^ ../includes/*.o

Make_config_mcc9_2D_muon_CC0pi: Make_config_mcc9_2D_muon_CC0pi.cpp 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(CC0PI_ANALYZER_FLAGS) -O3 -o $@ $^ ../includes/*.o

Make_Plots: Make_Plots.cpp 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(CC0PI_ANALYZER_FLAGS) -O3 -o $@ $^ ../includes/*.o
	
Make_Plots_2D: Make_Plots_2D.cpp 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(CC0PI_ANALYZER_FLAGS) -O3 -o $@ $^ ../includes/*.o
	
Make_cc0pi_effpur_pmu: Make_cc0pi_effpur_pmu.cpp 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(CC0PI_ANALYZER_FLAGS) -O3 -o $@ $^ ../includes/*.o

Make_cc0pi_effpur: Make_cc0pi_effpur.cpp 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(CC0PI_ANALYZER_FLAGS) -O3 -o $@ $^ ../includes/*.o
		
Make_Plots_test: Make_Plots_test.cpp 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(CC0PI_ANALYZER_FLAGS) -O3 -o $@ $^ ../includes/*.o

tutorial_slice_plots_test: tutorial_slice_plots_test.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^ ../includes/*.o
	
migration_matrix_cc0pi: migration_matrix_cc0pi.cpp
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(CC0PI_ANALYZER_FLAGS) -O3 -o $@ $^ ../includes/*.o	
	
.PHONY: clean

.INTERMEDIATE: stv_root_dict.o

clean:
	rm -f $(wildcard *.o) $(patsubst %.cpp, %, $(wildcard *.cpp))  Make_Plots Make_Plots_2D Make_Plots_2D_slicing_test 