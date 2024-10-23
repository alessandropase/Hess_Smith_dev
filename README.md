# Aerodynamic analysis of landing procedure using autoproducted Hess-Smith method 
## Description
This project implements the analysis of the landing procedure using an autoproducted Hess-Smith code. The purpose is to analyse the ground effect on a three surface wing while it approaches the ground.
The analysis is performed in 2D, considering the airflow as inviscid and irrotational, using Hess-Smith panel method implemented in MATLAB. The code has been validated comparing the results of the isolated foils with the results obtained through XFoil.

The advisor of the project was Dr. Andrea Rausa, teaching assistent of the Aerodynamic course of professor Franco Auteri.

## Requirements
- **MATLAB**
- **XFoil** for validation of numerical methods


## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/alessandropase/NonLin_Mech_analysis.git
   ```
2. Ensure you have MATLAB installed along with the required toolboxes.


## Project Structure
- `Data for panelization of airfoils` — Panel definitions of used airfoils obtained using XFoil;
- `main.m` — Main code where the analysis is performed;
- `Functions` — MATLAB Functions that are necessary for the analysis;
- `README.md` — Project documentation.

## External Resources
- [XFoil](https://web.mit.edu/drela/Public/web/xfoil/) — Software used for validation.

## Acknowledgment
This project was developed with Alessandro Pasolini.

  
