# ATLAS_2023_BJETS  
**Monte Carlo Simulation of High-Energy Particle Collisions**

*University of Glasgow Physics Honours Project - PHYS4022P*

---

## Project Overview
**Title:** *"Measuring b-Quark Jet Structure at the LHC"*  
**Supervisor:** Dr. Andy Buckley  

This project analyzes bottom quark (`b-jet`) structures in proton-proton collisions at the LHC, focusing on top quark pair (`tt̄`) decay channels. The study investigates jet substructure observables including:
- **N-Subjettiness**
- **Les Houches Angularity**
- **Energy Correlation Functions (ECFs)**
- **C₂ and D₂ correlations**

### Key Objectives
1. Develop a "Tag and Probe" methodology for `b-jet` identification
2. Quantify background contamination from light quarks (u, d, s) and gluon jets
3. Compare simulation results with ATLAS experimental data
4. Evaluate event selection efficiency and purity metrics

---

## Technical Implementation
### Core Components
- **Simulation Framework:** PYTHIA 8.3 + Rivet 3.1.6
- **Event Generation:**
  - **Signal Process:** `tt̄ → WWbb → ℓνℓνbb` (dileptonic decay)
  - **Background Processes:** 
    - WW/ZZ boson production
    - Z+jet events

### Dependencies
| Software       | Installation Guide                     |
|----------------|----------------------------------------|
| Docker         | [docker.com/get-started](https://www.docker.com/get-started) |
| Rivet          | [rivet.hepforge.org](https://rivet.hepforge.org/) |
| PYTHIA         | [pythia.org](https://pythia.org/)      |

---

## Workflow Setup
### 1. Docker Environment Configuration
```bash
# Start Rivet-PYTHIA container with directory mounting
docker run -it --rm -v $PWD:/host hepstore/rivet-pythia
# Compile analysis code
rivet-build ATLAS_2023_BJETS.cc -o RivetATLASBJETS.so
# Generate events (example: 10k tt̄ events)
pythia8-main93 -c ttbar-dilep.cmnd -o results.yoda -n 10000

