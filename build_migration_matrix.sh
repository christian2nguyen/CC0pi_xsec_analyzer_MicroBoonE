g++ -o migration_matrix migration_matrix.cpp \
  `root-config --libs --cflags` -lEG -lEGPythia6 -lGenVector -lGeom

g++ -o migration_matrix_cc0pi migration_matrix_cc0pi.cpp \
  `root-config --libs --cflags` -lEG -lEGPythia6 -lGenVector -lGeom

