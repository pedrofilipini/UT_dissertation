#######################
#multibart folder
#######################
Contains a version of multibart that can run the scalable tsBART/tsbcf variant.

#######################
#Auxiliary folder
#######################
Contains the gp_approx_fun.R file, which has functions that will later be added in the final R package.

#######################
#Dataset folder
#######################
Contains the file bcfstudy1_long.csv, which has the data that has been analyzed.


#######################
#study_1_model_2024.R
#######################
File used to "clean" the dataset and ran the model.
- Uses a scalable version of tsbcf (number of elements in t variable can be large);
- Length-scale is sampled within the model, just be careful setting the expected number of crossings if needed.



#######################
#study_1_graphs_2024.R
#######################
File used to generate the graphs of the report. Generates the graphs for:
- Curves over time (treatment and control);
- CATEs over time;
- CATEs by sex subgroup;
- CATEs by negative/positive prior mindsets (+ dual search for these subgroups).
A version of the model including self-esteem as moderator has been ran, but abandoned since the effects were low.


#######################
#tsbart_timing.R
#######################
Compares tsBART time with scalable tsBART.
Produces graph that shows the improved scalability with a toy example.

#######################
#ls_update_sim.R
#######################
Compares scalable tsBART with fixed and variable length-scales.
Produces graph that shows the improved results with a toy example.

#######################
#basi.R
#######################
Produces graphs in the appendix related to the Gaussian process approximation.

#######################
#error_graphs.R
#######################
Produces graphs in Section 2.4 related to minimizing the relative error.