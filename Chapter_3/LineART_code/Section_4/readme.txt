Description of the files/folders

Obs: Remember to change the directory before saving/running things.
Obs 2: I am using foreach for 20 threads, most computers cannot handle that. Check before running. 

------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

Dataset folder
-Contains 41 out of the 42 datasets from Chipman et al. (2010);
-x.txt files contains the covariates;
-y.txt files contains the response.

------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

sizes_paper.csv
-Contains details of the categorical variables for each dataset.

------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

Simulations_datasets_20241010_foreach_bart.R
Simulations_datasets_20241010_foreach_lbart.R
Simulations_datasets_20241010_foreach_softbart.R
-Originally a single file, separated in three due to the amount of time it took to run;
-The "bart" and "softbart" files run the default settings, while "lbart" has a considerable amount of options;
-This was done because a default for lineart was needed and that was a good way of exploring its behavior;
-Datasets are separated in training (5/6) and test (1/6)

------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

5fold.R
-5-fold cross validation performed on each of the available options;
-For all of those, the scale prior for mu was set to on;
-RMSE for training and test used the posterior mean;
-Used to have a CV version for each dataset.

------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

best_fold.R
-Used to find the setting with the lowest average RMSE over the 5-folds for each dataset.

------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

finding_default.R
-Used to analyze, among all setting, which settings won against bart the most;
-The winner is the new default.


------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

general_comparison.R
-Used to compare the defaults of bart, lineart, and softbart with lineart CV;
-The winner is the new default.


------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

comparison.R
-Has comparisons of lineart and bart defaults;
-Basically shows that, in general, the results are aligned, but when it is better, it is considerably better.


------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

timing.analysis.R
-Compare the time to run the defaults for softbart and bart in relation to lineart default.



