# Code for "Influence of near-surface currents on the global dispersal of marine microplastic"
Authors: David Wichmann, Philippe Delandmeter, Erik van Sebille 

This repository contains all code needed for the simulation and analysis for the paper. For questions, please contact d.wichmann@uu.nl.

# OceanParcels versions
For all but the 3D simulation, we used OceanParcels version 1.11 with slight modifications: 
For loading a sub-set of depths, it was necessary to comment lines 173 and 174 of field.py in the parcels folder.

We used OceanParcels version 2.0 for the 3D simulation. See: Delandmeter, P. & van Sebille, E. The Parcels v2.0 Lagrangian framework: new field interpolation schemes. Geosci. Model Dev. Discuss. 1–24 (2019).

# Particle Initial Coordinates
See InitialCoordinates/CreateGrid.py for creating a uniform 2D initial grid of particles.

# Simulation
For all but 3D simulation: Use Simulations/AdvectParticles.py
For 3D siulation: AdvectParticles_3D.py

## Simulation for fixed depths
The following command will put particles in the 13th depth level of the NEMO data. Initial particle locations are split into different sub-sets. The index posidx refers to this index. Here we use the initial grid no 3:

```python AdvectParticles.py -name FileName -y 2000 -m 1 -d 5 -simdays 3650 -posidx 3 -depth 13```

## Simulation for Uniform mixing
The following command will execute a simulation where particles are uniformly replaced along the vertical at each time step between 0 and 120 m.

```python AdvectParticles.py -name FileName -y 2000 -m 1 -d 5 -simdays 3650 -posidx 3 -uniformmixing True```

## Simulation for Kukulka 2012 mixing
Simulation for exponential mixing according to Kukulka 2012. The rise speed is set to 0.003 m/s.

```python AdvectParticles.py -name FileName -y 2000 -m 1 -d 5 -simdays 3650 -posidx 3 -kukulkamixing True -wrise 0.003```

## Simulation for 3D particles

python AdvectParticles_3D.py -name FileName -y 2000 -m 1 -d 5 -simdays 3650 -posidx 3

# Analysis
Execute Analysis/figures_paper.py to create all figures and tables for the main document and the supplementary material (adjust read-in file names according to the simulation output file names).  Analysis/AnaObjects.py contains the objects and region definitions for the data analysis.
