# TODO list 

# Figures 

Fig 1 
    - [] torch cross section with CD nozzle 
    - [] make jet flare
    - [] dimensions on schematic (jet size and goldilocks)
    - [] flip C, crop the top a bit. 

Fig 2

 - [x] CFD camera compare  and lifetime measurement
 - [] heat balance 

Fig 3 - 
    Comparison to CFD nK
    - [] add absorption measurement and fit

Figure 4
- [x] 2 positions 
- [x] remove C and PD


Figure 5 
- [] fix colors for position confusion
- [x] ionization process (only and B)

Figrue 6
recombination process 

- [x] add delta ne0


# Immediate

- [] noneq_modeling
    - [] implement unit handling noneq
    - [x] light absorption KOH case
    - [x] add eta based on experiment

    - [] double check new O2 kinetics
        - [] implement other species

 - [] why is oxygen not dominating Kp?


- [] get tau and ne_0 for nominal case. Calcaulate eta_experiment. Add to discussion

- [] fix CFD barrel. Think it has to do with spike. Check beam position. 

- [] move diagnostics to SI?

- [x] try fitting mws to a profile with O2 tau and K+ bimolecular combo
    - [] try refactoring dnedt models

- [x] lifetime comparison of 536 and 516. Need to compare across positions. Relevant location for 516 may be inward. 

- [] collect diagnostic information
 - [] oscilliscope settings time information
  - [] try taking off coarsening

- [x] investigate why phi coordinate value is 0.65 and understand motor coordinates


- [] add standard plots for all three sequences to SI
    - [] fit parameters (nK, tau, ne_0)
    - [] raw 