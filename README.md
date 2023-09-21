# This is a fork of main branch

# Cellular and Neurochemical Basis of Sleep Stages in the Thalamocortical Network

This project reproduces the results of [Cellular and neurochemical basis of sleep stages in the thalamocortical network](https://elifesciences.org/articles/18607). The study investigated how acetylcholine, histamine, GABA, and monoamines interact with shifting the brain between waking, non-rapid eye movement (NREM) sleep and rapid eye movement (REM) sleep within an ultradian cycle.

## Usage

The code is implemented in C++ and uses OpenMP multi-core parallelism to optimize the code's performance.
To generate the network connectivity file, navigate to the root directory of the project and run the following command: 

`make network network_config=network.cfg`

To run the simulation, use the following command:

`make run` 

The code generates several output files in the out folder that illustrate the electrical activity in the brain under different levels of neuromodulators. The output files include the membrane voltage of cortical neurons (time_cx) and other neuron types. You can vary the acetylcholine, histamine, GABA, and monoamine levels by modifying the parameters in the params.txt file. The network connectivity is given in the `network.cfg` file.

The `network.cfg` file includes the connectivity information of different brain regions and neuron types. The `params.txt` file contains several parameters, such as the levels of neuromodulators, the time step, and the simulation duration. You can modify these parameters to customize the simulation.

The output files are saved in the out folder with a timestamp in the filename. Each file contains the activity of different neurons at different times during the simulation. The membrane voltage of cortical neurons is saved in the time_cx file. Other neuron types are saved in separate files with their respective names.

## About the model

The computer model used in the study considers the known effects of each neuromodulator on different types of brain cells and how they interact with each other. The computer model successfully reproduces the patterns of brain electrical activity observed during different stages of sleep. The study found that acetylcholine, GABA, and monoamines work together to modulate the activity of different types of brain cells, leading to the transition between sleep stages. The model also predicted that the relative levels of neuromodulators could lead to characteristic brain EEG rhythms during different stages of sleep. The study also found that during NREM sleep, the power of spindle and delta oscillations is negatively correlated in humans and positively correlated in animal recordings. The differences in the relative level of acetylcholine explained this discrepancy.

The study provides insights into the complex interactions between neuromodulators and brain activity. It demonstrates the potential of computer models in investigating these interactions. The findings can contribute to developing new therapies for sleep disorders and other neurological conditions. The C++ code in this project can be used as a starting point for further research on the cellular and neurochemical basis of sleep stages in the thalamocortical network.

For more details on the study, please refer to the publication: Krishnan et al., Cellular and neurochemical basis of sleep stages in the thalamocortical network, eLife 5:e18607, https://doi.org/10.7554/eLife.18607.




