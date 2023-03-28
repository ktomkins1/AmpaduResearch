# Virginia Tech Optically Injected Feedback Laser System

Create figures for exploring the stability regions for laser parameters

To run, use `initiation.py`

The simulation is controlled by a Configuration File.  Configs are a set of options of the simulation.  A config file is a json representation   Additionally, a Model File is used as the mathematical basis for the simulation.

### Config Files
Configuration files are kept in `configs/` in the working directory.  For example, if run from the desktop, then `Desktop/configs/<config name>.json` should be defined.

Configuration files can be created in two ways:
- Edit `config_create.py` to change the values of the default configuration, making sure that the documentation stays in place, then run `config_create.py`
- Import `config_create.py` into a python script, or at the (i)python prompt, and modify the members of `config_create.config`, then call `config_create.save_config`

To specify a configuration file other than `configs/config.json`, use the `-c` command line argument to `initiation.py`.  In this case, you may reference a file in the `configs` directory, or you may specify the path to the file.  Paths containing "configs/" will be assumed to be a path to a file in `configs/`.

A copy of the config used will also be present within the results.

Note to maintainers:  As much as possible, it is preferred that the system is backwards compatible at least a few "features" back.  Therefore, an older config should be able to be used with no difficulty.  For example, if a new configuration parameter "my_new_param" is added to the config, and an older config is used, a KeyError would get thrown wherever the new feature was introduced.  If this is handled (i.e. using a try except block), an old config would not cause an error.

### Model Files
Model files contain the set of equations of the differential system being simulated.  These can vary slightly, such as for a change of basis or notation.  Alternatively, they could be entirely different equations.

Model files can be created for a system, so long as it represents a list of "functions".  The simulator uses SymPy to handle the differential set, and requires that the model supports a setup method which takes 5 real values as input and returns a list of functions.  Nominally, these 5 real values would be those of P, Delta, alpha, eta, and T.  This setup method would then return the specific set of diffEqs for that configuration of variables.

A model file is not required to return 3 functions, although that is the number which the system was designed for.  It also is known to work with 2-function models.

The model file is specified in the Config.  Model files are kept in the repository (not the working directory) in `models/`

### Theory of Operation

The simulation framework relies on SymPy for the symbolification of the diffEq functions, and SciPy's `ode` and `solve_ivp` libraries to run simulation steps.

The simulation handles automatically finding the bounds of simulation and utilizing multiple threads (if supported).

Essentially, the workflow of the framework is as follows:
- Prepare a configuration with a certain set of variables, with one of them sweeping
- Break the configuration into a set of static configs
- Dispatch these configs to threads, limited by the number of threads of the system
- The trace which is run uses these parameters, along with the functions from the model
- The results of the trace in the form of raw data gets packaged into a dictionary
- The results are also used to create a plot, using parameters from the config
- The plots, the associated results, and the associated config file are saved to disk in `results\`

### Command line options

The `initiation.py` module can be run from the command line or in a batch script and has the following options:
- `-c`: use a config file other than `configs/config.json`
- `-r`: replay a results file, unzipping the results and plotting them again
- `-z`: specify an alternative name for the results file

### Future work and How to Contribute

The framework could be improved in several ways, such as more fully featured plotting, post-simulation parsing, documentation.

There are several `TODO`s in the code, as well as some partially implemented features which don't work correctly and are disabled.

The addition of new models with alternate structures or parameters could necessecitate a paradigm shift in the simulation pipeline.

Originally, this framework was created by Kevin Tomkins, and he can be reached at tomkinskevin1@gmail.com
