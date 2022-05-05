## MetaPredict workflow

### General description of this project

This is an annotated workflow for training machine learning models to predict the presence or absence of [KEGG metabolic modules](https://www.genome.jp/kegg/module.html) in bacterial genomes recovered from sequencing efforts.

Genome datasets used for models in this project are from the NCBI RefSeq bacterial genomes database (version 209, 2021) categorized as "complete genomes", and the set of GTDB bacterial genomes (version r95, 2021) that were assessed to have completeness of 100, contamination of 0, strain heterogeneity of 0, and MIMAG quality of "High Quality" (and are not also present in the RefSeq genomes used). A total of 30,646 genomes are utilized here for training machine learning models and assessing their performance. 424 models have been trained and tested; one model has been created for each KEGG module which is present in at least 0.1% of these 30,646 genomes. All model training and assessment was done using [Tidymodels](https://www.tidymodels.org/).

The [rscripts/](rscripts) folder contains code scripts in the following categories: [training_scripts](rscripts/training_scripts), [model_assessment_scripts](rscripts/model_assessment_scripts), and [sql_db_scripts](rscripts/sql_db_scripts). [training_scripts](rscripts/training_scripts) has code for training machine learning models utilized in this project; [model_assessment_scripts](rscripts/model_assessment_scripts) contains scripts for assessing the quality of trained models, and how they would perform on new, unseen data; [sql_db_scripts](rscripts/sql_db_scripts) has code used to create a SQL database containing all trained models that MetaPredict will query to make predictions on user input data.

The [bash_scripts](bash_scripts/) directory has scripts for training and assessing machine learning models on WHOI's Poseidon HPC.

More details will be added here to give background information about KEGG metabolic modules, as well as the rationale for this project and the model types used. Plots and descriptions of performance metrics for MetaPredict's models will also be posted here.
<br>  
<br>  

## Details about the R scripts in the rscripts/ directory

### Code in the [rscripts/training_scripts](rscripts/training_scripts) directory:
  - [create_model_env_vars.R](rscripts/training_scripts/create_model_env_vars.R): this script contains code to create train and test dataset splits for each model. Cross validation folds for nested cross validation are also created with this script; this is a list of various random samples of the training data that will be used to tune model hyperparameters. These training objects are then saved in a compressed binary format (.rdata) to be loaded in later model training scripts. If modeling is re-done later, this step is already complete and future modeling jobs can start with this pre-generated information.
 
  - [create_nnet_tune.R](rscripts/training_scripts/create_nnet_tune.R): code to assess a grid of possible hyperparameter values of a fully connected, multilayer perceptron neural network model on training data. The script saves the results of a grid search of the hyperparameter space via cross validation, emphasizing the following metrics: Cohen's Kappa and Receiving Operating Characteristic (ROC_AUC).
 
  - [create_xgboost_tune.R](rscripts/training_scripts/create_xgboost_tune.R): code to assess a grid of possible hyperparameter values of an Extreme Gradient Boosting (XGBoost) model on training data. The script saves the results of a grid search of the hyperparameter space via cross validation, emphasizing the following metrics: Cohen's Kappa and ROC_AUC.
 
  - [train_ensemble.R](rscripts/training_scripts/train_ensemble.R): code that fits a LASSO model to the combinations of hyperparameter values from each cross validation grid search iteration of **create_nnet_tune.R** and **create_xgboost_tune.R** jobs. This LASSO model will incorporate or eliminate models into a model ensemble based on the performance of the hyperparameter values that were assessed during cross validation. Models incorporated into an ensemble will be weighted according to how well they did on the hold-out folds during cross validation. This will ultimately create a group of models that will make predictions (in terms of probability) about the presence or absence of a KEGG module in a bacterial genome. Each probability calculated by a model will then be multiplied by the coefficient it was assigned by the LASSO model. The final probability computed by an ensemble model is a blend of all the weighted probabilities calculated by the model's ensemble members.
 
  - [fit_xgboost_newest_april_2022.R](rscripts/training_scripts/fit_xgboost_newest_april_2022.R): this code will fit an individual XGBoost model, using the hyperparameters from **create_xgboost_tune.R** that had the highest Cohen's Kappa and ROC_AUC scores. Performance of these models on held-out test data can then be compared to the performance of the ensemble models trained using **train_ensemble.R**.
 
  - [ensemble_master_script.R](rscripts/training_scripts/ensemble_master_script.R): this is a compilation of all the training scripts menteioned above that are used to train a neural network/XGBoost ensemble model.
<br>  

### Code in the [rscripts/model_assessment_scripts](rscripts/model_assessment_scripts) directory:
- [ensemble_test_on_simulated_fragmented_genomes.R](rscripts/model_assessment_scripts/ensemble_test_on_simulated_fragmented_genomes.R): This code assesses ensemble model performance (via model metrics including ROC_AUC, F1, Precision, Recall, Positive Preciction Value, Negative Prediciton Value, etc.) on simulated fragmented genomes. The 50 genomes used for testing here were assembled from quality-trimmed sequencing reads that were randomly downsampled in increments (10% - 90% of reads retained, in 10% increments; also 0.1%, 0.4%, 0.7%, 1%, 3%, 5%, 7% reads retained). These 50 genomes are a second held-out test dataset not included in model training and testing that are specifically for this purpose of testing models on simulated genomes of which only fragments have been recovered. The script saves the results of these model assessments.

- [xgboost_test_on_simulated_fragmented_genomes.R](rscripts/model_assessment_scripts/xgboost_test_on_simulated_fragmented_genomes.R): the code has the same purpose as **ensemble_test_on_simulated_fragmented_genomes.R** but is tested on XGBoost models, and the model assessment results are then saved.

- [compare_xgboost_to_ensemble_model_metrics.R](rscripts/model_assessment_scripts/compare_xgboost_to_ensemble_model_metrics.R): code to compare ensemble and XGBoost model performance on simulated fragmented genomes. Creates a multifaceted plot that shows differences in various model performance  metrics between the two model types.

- [assess_which_ensemble_jobs_need_to_run_on_poseidon.R](rscripts/model_assessment_scripts/assess_which_ensemble_jobs_need_to_run_on_poseidon.R): based on job logs on the Poseidon HPC, this code is simply determining which jobs failed and need to be re-run, or still have not been run yet. 
<br>  

### Code in the [rscripts/sql_db_scripts](rscripts/sql_db_scripts) directory:
- [create_metapredict_model_sql_database.R](rscripts/sql_db_scripts/create_metapredict_model_sql_database.R): In the MetaPredict package, the models will be stored in an SQL database and will either be sequentially or batch loaded into memory to make predictions on user input data, or predictions will be run inside the database without loading any of the models into memory. The top sections of the code contain various scratch code used to create a sample SQL database and test querying the test database. The final bottom code section titled "create official MetaPredict SQL database" is where the current SQL database containing all current ensemble models was constructed on Poseidon.
<br>  
<br>  

## Details about the bash scripts in the bash_scripts directory:

### Code in the [bash_scripts](bash_scripts/) directory:
- [compute_train_all_ensemble_array.sh](bash_scripts/compute_train_all_ensemble_array.sh): an array slurm script that calls training R scripts to train ensemble models, one per compute node. Sequentially calls create_model_env_vars.R, create_nnet_tune.R, create_xgboost_tune.R, and train_ensemble.R. Output is saved from each script so that future iterations of model training can be started up with intermediate output data from any of the scripts. 

- [scav_train_all_ensemble_array.sh](bash_scripts/scav_train_all_ensemble_array.sh): same as compute_train_all_ensemble_array.sh but it's executed on scavenger nodes (instead of regular compute nodes).

- [scav_only_fit_ensembles_array.sh](bash_scripts/scav_only_fit_ensembles_array.sh): code that only calls train_ensemble.R, and uses pre-generated .rdata files created using create_nnet_tune.R and create_xgboost_tune.R.
