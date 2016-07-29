# FocalNet - Estimation of focal lengths from correspondences #

## Important documents ##
 - **Learning algebraic geometry - notes**: Main document with the ideas and thoughts documented, keeps track of the research. https://docs.google.com/document/d/12KifffpP7B1aR7QZmrOxUjuLTH_8IvQM4KGmNIk-n5M/edit?usp=sharing
 - **Learning algebraic geometry - results**: Document with partial results containing a lot of graphs and data. https://docs.google.com/a/fel.cvut.cz/document/d/1SvP7J-ywYK8i5ygudt7YVuBnE6yvyaVFwagj8vLB8Mo/edit?usp=sharing
 - **LaTeX document**: Document describing main ideas. https://cs.sharelatex.com/project/577cc312485e058373049940

## Git repository ##
To keep track of the changes in the code, the Git tool is used. Data folders are excluded from this repository because of the large size of the files. Please do not forget to set up your name and e-mail before first commit: `git config --global user.name "John Doe"`, `git config --global user.email johndoe@example.com`.

## Folder structure ##
 - `deepag/`
   - `bin/`
     - `deepag\`: Main folder with the most important scripts.
       - `results\`: Folder containing results saved by the scripts. The results can be viewed by the `plot*` scripts.
       - `deepag.m`, `deepag_parallel.m`: Scripts used to generate the training, validating and testing datasets from the raw data. The results are saved in `data/` folder as `features.mat`, which contains generated feature vectors and as `correspondences.mat`, which contains selected correspondences.
       - `cons_predictor.m`, `closed_form.m`: Scripts to compute error of the constant predictor and of the closed form linear regressor.
       - `k_nn.m`, `k_nn_GPU.m`: Scripts to compute the error of the k nearest neighbours regressor in the feature vector space. The version `_GPU` is designed to use on the GPU. 
       - `k_nn_pixels.m`, `k_nn_pixels_generateData.m`: Scripts to compute the error of the k nearest neighbours regressor in the image coordinates space. `_generateData` version is used to generate synthetic data for testing. The generated data are saved in `data/` folder as `correspondences_syntetic.mat`.
       - `nn*`: Scripts with neural network regressors. The scipts differ in the number of non-linear layers. The `_GPU` versions are used to be executed on the GPU.
       - `plotKNN.m`, `plotNN.m`: Scripts to plot and inspect the results from `k_nn*` and `nn*` scripts respectively. Very handy when the scripts are executed on the cluster on which you can not plot anything.
     - `lib`: Folder with libraries, like `MatonvNet` and utils functions.
   - `data/`: Folder containing all data files.
 - `export/`: Folder with exported images and other documents, usually used in the gDoc with partial results.
 - `logs/`: Folder with log files from batch jobs executed on the cluster
 - `coef_mat_two_focal.m`: Matlab function generating feature vector from seven correspondences as Zuzana sent it.
 - `*.pbs`: Script filed used to submit jobs to the cluster.
 - `Whiteboard-*.jpg`: Captured whiteboards after some discussions.
