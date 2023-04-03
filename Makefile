all: mindcurrent

mindcurrent: main.cpp CellSyn.cpp CellSyn.h currents.cpp currents.h io.cpp io.h network.h network.cpp
	g++ -g -O3 -lm -Wall -fopenmp  main.cpp CellSyn.cpp io.cpp currents.cpp network.cpp -o mindcurrent
	 #icpc -O3 -Wall -fopenmp main.cpp CellSyn.cpp io.cpp currents.cpp network.cpp -o mindcurrent

#sanity check whether current version of code changes output (compared to previously stored test/.files)
check: mindcurrent
	./mindcurrent test/params.txt test connection_info2
	cd test; for f in *; do diff -u $$f .$$f; done

check-prepare: mindcurrent
	./mindcurrent test/params.txt test connection_info2
	cd test; for f in *; do cp -f  $$f .$$f; done

run: mindcurrent
	./mindcurrent params.txt out connection_info2

doxy:
	doxygen ./docs/Doxyfile

clean: 
	-rm mindcurrent 

network:
	g++ -O2 generate_network.cpp -o generate_network
	./generate_network $(network_config) $(mri_network) $(3D_subnet) $(3D_distance)> connection_info2



  