# TODO list 

# Immediate

- [] noneq_modeling
    - [] implement unit handling noneq

    - [] double check new O2 kinetics
        - [x] implement other species


- [] get tau and ne_0 for nominal case. Calcaulate eta_experiment. Add to discussion

- [] fix CFD barrel. Think it has to do with spike. Check beam position. 

- [x] collect diagnostic information
 - [] oscilliscope settings time information
  - [x] try taking off coarsening

- [] add standard plots for all three sequences to SI
    - [] fit parameters (nK, tau, ne_0)
    - [] raw 

- [x] Input properties table. 

- [x] add phi = 0.6 to simulation input report. 

- [] add calorific value and CH ratio to sim report
 - [] reprocess 


 change particle to molecule


- [] btw, the case window for 0.99_1 case on 5/12 could be shortened. The fuel is still changing at the start of the window. If you start window at  11:29:55(10 seconds later)  it will be better. though I'm just looking at DAQmx data right now. The fuel pump and combustion pressure change behavior half way through this case, so I wonder if something with syringe pump changes.Change 

- [] syringe zero removal duplicated in post processing? 

- [] 3% case on 05-12

- [] make pre norm cutoff None by default (And similar kwargs) and update/test pipelines