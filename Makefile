bin/bayesianfastalign: src/bayesianfastalign.cc
	g++ -std=c++11 -Wall -g -O3 -L/home/austinma/prefix/lib -I/home/austinma/git/cpyp/ src/bayesianfastalign.cc -o bin/bayesianfastalign -lboost_program_options -lboost_serialization
