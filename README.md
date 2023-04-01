# Sleep Stages  

Code to replicate results from the paper : 
Cellular and neurochemical basis of sleep stages in the thalamocortical network. eLife 5:e18607. https://doi.org/10.7554/eLife.18607

The code is written in C++ and uses openmp for running in multiple cores. 

To run the simulation from shell:
```sh
#Create output folder
mkdir out

#Generate network connectivity file
make network

#Generate network connectivity file
make network

#Run the simulation
make run
```

There will several files generated in the output folder (/out) which includes the membrane voltage of cortical neurons (time_cx) and other neuron types. 

Network connectivity is given in network.cfg. Many other parameters are specified in params.txt
