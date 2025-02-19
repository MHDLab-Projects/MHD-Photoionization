# Final Data and Figures

Directory for final dataset and figures in paper. The intention is cleaned up scripts outputting datasets that are then used to create select image files without exploratory analysis, like in `experiment/analysis` folder. Note that figures used the SI may still be produced in other folders.  

* dataset - scripts to generate final dataset including analysis (fitting/processing/etc.). These steps are often time-intensive so this allows for iterative development of figures afterward. 
* figure_panels - scripts to produce plots from final dataset. Sometimes these scripts output fully laid out figures into the figures directory, or individual panels that will be laid out in manually created svg files. 
* figures - svg files that layout figure panels into final figures
* analysis - additional analysis of final dataset for supplementary figures beyond paper figures.