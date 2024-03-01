# TODO list 

# Figures 

Fig 1 
    - [] torch cross section with CD nozzle 
    - [] make jet flare
    - [] dimensions on schematic (jet size and goldilocks)
    - [] flip C, crop the top a bit. 

Fig 2

 - [] CFD camera compare  and lifetime measurement
 - [] heat balance 

Fig 3 - 
    Comparison to CFD nK
    - [] add absorption measurement and fit

Figure 4
- [] 2 positions 
- [] remove C and PD


Figure 5 
- [] fix colors for position confusion
- [] ionization process (only and B)

Figrue 6
recombination process 

- n_e fit curve
taur and delta n_e

Figure 7
 = [] why is oxygen not dominating Kp?

# Immediate

- [x] respond comments
- [] move diagnostics to SI?
- [] update noneq in SI
- [] run through writing
- [x] add summary of KOH model (email)

- [x] fix MWS data based on no torch transmission. What dates do we have No Torch MWS data for? can we check it is the same as far downstream. 
    - [] calculate additional statistics of mws position data. Probably don't need T0 to explain decrease in AS. Compare using T0 to pre pulse magnitude average. 
# Medium

- [x] implement numerical integration nK determination to final pipeline
- [] improve main pimax image
    - [x] calibration
    - [x] compare width to cfd beam profile, also compare shock diamonds
    - [ ] double check dimensions
    - [ ] reduce downsampling. 

- [x] make AS plot with n_e as y axis
    - [] all dates for main plot

- [x] use beam profile, CFD KOH, and KOH absorption cross seciton to estimate delta n_e

- [x] laser profile. 


- [x] try fitting mws to a profile with O2 tau and K+ bimolecular combo
    - [] try refactoring dnedt models
- [] double check new O2 kinetics
    - [] implement unit handling noneq
    - [] implement other species

- [] add phi=0.6 and KOH data to mws position plot. 

- [] lifetime comparison of 536 and 516. Need to compare across positions. Relevant location for 516 may be inward. 

- [] collect diagnostic information

# Long

- [] investigate why phi coordinate value is 0.65. 
