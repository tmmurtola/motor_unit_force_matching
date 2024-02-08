# motor_unit_force_matching

This package contains MATLAB code that generates motor unit (MU) pools and uses them to perform force matching tasks used in

Murtola & Richards (2024) Matching dynamically varying forces with multi-motor-unit muscle models: A simulation study [submitted].

# Preliminaries

Requirements:
- Matlab (2023a), compatibility with other recent versions is likely but not guaranteed
- Parameter optimisation also requires Matlab toolboxes: Global optimisation toolbox (4.8.1) and (optionally) Parallel computing toolbox (7.8)

Package contents:
- run_force_matching.mlx - main MATLAB script for generating MU pools and running force matching tasks
- optimise_neural_drive_parameters.mlx - MATLAB script to optimise and save neural drive parameters
- plot_firing_rate_figures.mlx - MATLAB script to simulate trapezoidal tasks and plot firing rate figures
- helper functions and class definitions:
	- MUpool - class definition for MUpool objects
	- muscle_MU_settings - function setting MU pool parameters
	- generate_MUpool - function generating a set of MU pools from parameters
	- force_matching_tasks - function generating force and FVL profiles
- reaching_data.mat - data for simulated reaching tasks
- optim_params_TA_precomp.mat & optim_params_SF_precomp.mat - precomputed set of neural drive parameters


# Usage

Make sure all functions and data files are on Matlab's search path. 

To use run_force_matching.mlx:
1. Use commands on lines 1-6 to select muscles and corresponding parameter files.
2. (optional) To modify MU pool parameters, open muscle_MU_settings.m and edit values on lines 8-41.
3. Run the script.


To use optimise_neural_drive_parameters.mlx:
1. Set parameters for the set of elementary tasks to be used for optimisation. 
	taskDur = duration of muscle length change (e.g. 2 sec), 
	initLength = the initial normalised muscle length (e.g. 1 for task to start at optimal length, l0)
	lengthChange = relative length change for each task (e.g. [-.2 0 2] for three tasks: concentric to 0.8*l0, isometric, and eccentric to 1.2*l0)
2. (optional) Change muscle-specific force scaling for tasks on lines 15/17, if desired.
3. Select muscle on line 7. (optimise one muscle at a time).
4. (optional) The optimisation uses precomputed parameter values as a partial initial population. To change the file used, edit line 35. To optimise without 
   initial population, comment out lines 35-36 and 49, and remove the variable initVals from the function call on line 55.
5. Set filename for saving results on line 62.
6. Run the script.



To use plot_firing_rate_figures.mlx:
1. Use first two sections similarly to run_force_matching.mlx.
2. For the third section, download reference data as instructed and specify its loaction.
3. Run the script.



