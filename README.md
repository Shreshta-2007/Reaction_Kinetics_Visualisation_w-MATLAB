# Reaction_Kinetics_Visualisation_w-MATLAB
# Reaction Kinetics Visualisation (MATLAB)

## Problem statement

How fast does a reaction really go, and how much does a few degrees change that?

This project turns Decomposition of dinitrogen peroxide's published kinetic data into plots you can read at a glance. Starting from first-order decay, C(t) = C0·exp(−kt), it shows how concentration falls for different rate constants. It then uses the Arrhenius equation to show how temperature changes k and the half-life. [A short thermodynamics section covers how the equilibrium shifts with temperature.]

##Reaction
N₂O₅(g) → 2 NO₂(g) + ½ O₂(g) 

Every parameter comes from literature, not made-up numbers, so the trends can be compared against real behaviour.

## Full report

for full analysis and all details  see [Reaction Kinetics MATLAB project report.pdf](https://github.com/user-attachments/files/34208359/Reaction.Kinetics.MATLAB.project.report.pdf)


 
## What the code does

- Plots concentration against time using C(t) = C0·exp(−kt) for several rate constants
- Uses the **Arrhenius equation** to calculate k at different temperatures
- Plots how k and the half-life change with temperature
- [Equilibrium part: for example, plots how the equilibrium constant changes with temperature using the van 't Hoff equation]
- All plots are labelled and saved to `figures/`

## How to run

1. Open MATLAB and set the project folder as the current folder.
2. Run `[main_script_name].m`.

Tested on MATLAB [version]. No toolboxes needed. [Delete this line if that isn't true.]

## Results

<img width="913" height="476" alt="Screenshot 2026-10-30 at 8 41 06 PM" src="https://github.com/user-attachments/assets/523d7651-f3e8-4e99-ae6c-3b69fab7237f" />
<img width="980" height="477" alt="Screenshot 2026-10-30 at 8 40 50 PM" src="https://github.com/user-attachments/assets/78cd2c69-aba5-4eaa-99d3-dd958a552c96" />



## Assumptions

- First-order, [irreversible/reversible] reaction
- Constant temperature for each run


## Author

Shreshta Kanduri  | BITS Pilani Dubai B.E. Chemical Engineering  | www.linkedin.com/in/shreshta-kanduri-311230312
