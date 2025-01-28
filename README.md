# Classical-Data-in-Quantum-Machine-Learning
This repository contains the code for the paper 'Classical Data in Quantum Machine Learning Algorithms: Amplitude encoding and the Relation between Entropy and Linguistic Ambiguity'.

This repository is structured as follows. 
There are three folders, corresponding to the three contributing sections in the paper. 
In the 'datasets' folder, the datasets (both new and original) are collected. 

In folder 'Replication', the code for the replication task is collected. The hyperparameters in the training procedure are set to one specific combination for all contributions. 

In folder 'AmplitudeEncoding', the code for the task of amplitude encoding classical data on a quantum circuit is collected. Loss curves, as well as various performance metrics are created after the training. 

Lastly, in folder 'Entropy', the code capturing the entropy of different density matrices modelling sentences with varying levels of ambiguity are presented. 

In 'modelParameters', the parameters of the trained models are saved as .pkl files. 
In 'weights', the weights of the different models are captured to be able to construct the VQCs after the training procedure, using Python. 
In 'blochPlots', the plots of the noun states are displayed on the Bloch sphere for the different models. 
