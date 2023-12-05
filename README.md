# ATLAS_2023_BJETS
Monte-Carlo Simulation of high-energy paticle collision

## University of Glasgow Graduation Project (PHYS4022P)
### "Measureing b-quark structure at the LHC", Supervised by Dr. Andy Buckley
This project is dedicated to measure in-depth the bottom quark jet structures, where jets are defined as collimated sprays of hadrons produced from high-energy collisions. The analysis includes various key characteristics of these jets, including the N-Subjettiness, Les Houches Angularity, Energy Correlation Functions, C2 and D2 correlations and Lund Jet Planes. The primary source of the bottom quarks in this project will be coming from the decay pattern of the top and antitop quarks pair (ttbar). This is chosen specifically due to its distinct decay channel of b-jets plus dileptonic, which serves as a complex enough yet distinguishable experiment.

High energy collisions yield a significant quantity of light quarks (u, d, s quarks) and gluon jets, known as the background contamination. A key emphasis in this project lies in establishing solid selection criteria for identifying valid b-jets: A method often referred to as the "Tag and Probe" method. Through the evaluation of the efficiency and purity of these events, the objective is to compare with real data from ATLAS at LHC for the detection of b-quarks.

## Implementation
The simulation is run via the Rivet and Pythia systems, respectively. Installation and general information can be found on their official website:
- https://rivet.hepforge.org/
- https://gitlab.com/hepcedar/rivet
- https://pythia.org/latest-manual/Welcome.html

Rivet is commonly run on Docker desktop, the installation of it can be found at its official page:
- https://www.docker.com/
<br />
Within the code files, comments and references are listed for easy understanding and modification. There are two major events simulation in the codes:
- Bottom quark signal: $ttbar$
- Background contaminations: $WW$, $ZZ$, $Z+jet$

In order to switch between the modes, modify ttbar-dilep.cmnd file by commenting out or turn on/off each condition. The detailed applications of the functions and code are all in the pythia manual. The ATLAS_2023_BJETS.cc file can be modified at own convinience for investigating different structures and filtering conditions.

## Startup
The procedure used to start-up a simulation is as follows:
1. On your computer with docker desktop installed and opened, insert
     - "docker run -it --rm -v $PWD:/host hepstore/rivet-pythia"
   in CMD, where $PWD is the directory where you would like to store the codes and simulation. You should be in the docker workplace if it's successfully pulled.
3. Having both the ATLAS-2023-BJETS.cc and ttbar-dilep.cmnd in the previous directory, type in
     - "rivet-build ATLAS-2023-BJETS"
   in the workplace and a shared object called RivetAnalysis.so (defult name, can be assigned) will be created. 
4. Start the simulation by typing
     - "pythia8-main93 -c ttbar-dilep.cmnd -o FILENAME -n EVENTNUM"
   where "FILENAME" will be the name of the output .yoda file, and "EVENTNUM" is the number of events to be run.

