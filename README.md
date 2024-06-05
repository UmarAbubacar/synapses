# Pyramidal Cell + Synapses Demo

Note: To use this demo you must be on the ```umar_neurite_test``` branch of BioDynaMo.

This demo simulates the growth of a pyramidal cell. The simulation consists of the following files:

The files in the `src` directory contain the implementation of the simulation.
synapses.h and synapses.cc are the files that contain the implementation of the simulation with the code for the custom neurons and synapses in the file basic_neuron. 
The synapse operation in this example is a behaviour defined in the synapse_op.h file that is scheduled to run on the last iteration of the simulation.


To compile and run the simulation, execute the following command in the terminal.

```
bdm run
```

The simulation will use the `bdm.json` parameter file and will create
visualization files in directory `output/synapses`.

To render an image of the final neuron, execute:

```
bdm view
```
Press play and interactivly adjust the camera to get the desired view of the neurons. 

```
paraview output/synapses/synapses.pvsm
```

For more information about the simulation itself have a look at the following publication:

Lukas Breitwieser et al. BioDynaMo: a modular platform for high-performance agent-based simulation.
Bioinformatics, 2021. DOI: https://doi.org/10.1093/bioinformatics/btab649.
