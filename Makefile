# Variables
CXX = g++
CXXFLAGS = -g -O3 -Wall -fopenmp -lm
TARGET = mindcurrent
SOURCES = main.cpp CellSyn.cpp io.cpp currents.cpp network.cpp
HEADERS = CellSyn.h currents.h io.h network.h

# Default target
.PHONY: all
all: $(TARGET)

# Build the main target
$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

# Check correctness against stored test files
.PHONY: check
check: $(TARGET)
	./$(TARGET) test/params.txt test connection_info2
	cd test; for f in *; do diff -u $$f .$$f; done

# Prepare baseline test files
.PHONY: check-prepare
check-prepare: $(TARGET)
	./$(TARGET) test/params.txt test connection_info2
	cd test; for f in *; do cp -f $$f .$$f; done

# Run the main executable
.PHONY: run
run: $(TARGET)
	./$(TARGET) params.txt out connection_info2

# Generate documentation
.PHONY: doxy
doxy:
	doxygen ./docs/Doxyfile

# Clean up build files
.PHONY: clean
clean:
	-rm -f $(TARGET) generate_network

# Build the network binary
.PHONY: network
network:
	$(CXX) -O2 generate_network.cpp -o generate_network
	# Uncomment the following line to run the network binary during the build
	# ./generate_network $(network_config) $(mri_network) $(3D_subnet) $(3D_distance) > connection_info2