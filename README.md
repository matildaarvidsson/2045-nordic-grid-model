# PyPSA-Eur
## Installation
This is an adapted version of the 'installation' page of the [PyPSA-Eur documentation](https://pypsa-eur.readthedocs.io/en/latest/installation.html)

The subsequently described installation steps are demonstrated as shell
commands, where the path before the `%` sign denotes the directory in
which the commands following the `%` should be entered.

### Clone the Repository
First of all, clone the adapted [PyPSA-Eur
repository](https://github.com/niekbr/pypsa-eur) using the version
control system `git`. You can either use the command line tool [Git Bash](https://git-scm.com/download/win) or a GUI like [GitHub Desktop](https://desktop.github.com).
This installation guide only shows the command line commands, but there is an alternative.

``` bash
cd /some/path
git clone https://github.com/niekbr/pypsa-eur.git
```

Checkout the `dnv` branch.
``` bash
cd pypsa-eur
git checkout dnv
```
### Install Python Dependencies

PyPSA-Eur relies on a set of other Python packages to function. We
recommend using the package manager and environment management system
`conda` to install them. Install
[miniconda](https://docs.conda.io/en/latest/miniconda.html), which is a
mini version of [Anaconda](https://www.anaconda.com/) that includes only
`conda` and its dependencies or make sure `conda` is already installed
on your system. For instructions for your operating system follow the
`conda` [installation
guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

The python package requirements are curated in the
[envs/environment.yaml](https://github.com/PyPSA/pypsa-eur/blob/master/envs/environment.yaml)
file. The environment can be installed and activated using

``` bash
conda env create -f envs/environment.yaml
conda activate pypsa-eur
```

Note that activation is local to the currently open shell! After opening
a new terminal window, one needs to reissue the second command!

> **Note**
> If you have troubles with a slow `conda` installation, we recommend to install [mamba](https://github.com/QuantStack/mamba) as a fast drop-in replacement via
> ```
> conda install -c conda-forge mamba
> ```
>
> and then install the environment with
> 
> ``` bash
> mamba env create -f envs/environment.yaml
> ```

###  Set Up the  Configuration

PyPSA-Eur has several configuration options that must be specified in a
`config.yaml` file located in the root directory. A good configuration for the DNV Nordics model is already included in
the configuration file `config.dnv.yaml`. It might be necessary to update some config based on your needs. 
More details on the configuration options are in the [PyPSA-Eur documentation](https://pypsa-eur.readthedocs.io/en/latest/configuration.html)

Before first use, create a `config.yaml` by copying the DNV config example.

``` bash
cp config.dnv.yaml config.yaml
```

## Usage
PyPSA-Eur uses the workflow tool `snakemake`. 

PyPSA-Eur suggests simplifying the network (aggregration of buses) to make solving easier, but this project requires the full-resolution model. It can be obtained as follows.
``` bash
snakemake -c all networks/elec.nc
```
This tells snakemake to run the required scripts to obtain the file `elec.nc`, which is the full-scale model in their own format. This format will serve as input for the DNV Nordics model as described below.

# Nordics model
## Installation
It is assumed that `git` and `conda` are already installed, as shown in the PyPSA-Eur installation instructions above.

### Clone the Repository
> **Note**
> Instead of cloning the GitHub repository, it is also possible to open the model in the 'Grid models based on open sources - Documents' Team.

As this is a private GitHub repository, request access to the GitHub first.
``` bash
cd /some/path
git clone https://github.com/niekbr/dnv-nordics-model.git
```

### Install Python Dependencies
Similar to PyPSA-Eur, create a conda environment, which automatically installs dependencies.

``` bash
conda env create -f envs/environment.yaml
conda activate dnv-nordics-model
```

### Integration with PSS/E
It is assumed that PSS/E is already installed.

> **Note**
> You need to have the licence USB inserted when installing PSS/E.
>
> If the program crashes on startup, uninstall it and reinstall with the licence USB inserted.

The python integration of PSS/E `psspy` has to be accessible for python. This can be done by adding some paths to the [environment variables](https://helpdeskgeek.com/windows-10/add-windows-path-environment-variable/). Ask colleagues for help if required.

PATH:
- C:\Program Files\PTI\PSSE35\35.4\PSSPY39
- C:\Program Files\PTI\PSSE35\35.4
- C:\Program Files\PTI\PSSE35\35.4\PSSBIN
- C:\Program Files\PTI\PSSE35\35.4\PSSLIB

PYTHONPATH:
- C:\Program Files\PTI\PSSE35\35.4\PSSPY39
- C:\Program Files\PTI\PSSE35\35.4\PSSBIN
- C:\Program Files\PTI\PSSE35\35.4\PSSLIB

Lastly, in order to use the bus coordinates when drawing the slider, check the setting in `Edit` > `Preferences` > `Diagram` > Use available geophysical information when auto-placing items

## Usage
Similar to PyPSA-Eur, the `snakemake` workflow tool is used. First, make sure the `elec.nc` file is present in the `input` directory. This file is the output of running PyPSA-Eur.

Check the configuration file `config.yaml` and update if required.

Then to build a PSS/E network for all snapshots with validation, run:

``` bash
snakemake -c all
```

Some tips for using `snakemake`
- The `snakemake` command always requires the `-c` flag with the amount of parallel cores it can use. By default, use `all` to run it fastest.
- If you do not want to run the whole workflow (e.g. to test something), for instance to only build the database, write the output file you want: `snakemake -c all database/2020-high/nordics.sqlite`. This stops the workflow after building the database.
- Snakemake will only run a certain script if its input files are newer than the output files. This speeds it up after the first run

### Losses
Since the load data includes losses, the (estimated) losses need to be subtracted from the load per bidding zone.
To get a good approximation of the losses, the configuration needs to be set by making an initial guess for the losses, building the PSS/E model, and solving the load flow to get the actual losses.
By updating the guessed losses with the resulting losses iteratively, it will converge to the right values in about 2/3 iterations.

The losses can be retrieved per bidding zone by opening the PSS/E file and creating a 'Area/Zone/Owner Total Report' ![aoz-icon.png](aoz-icon.png) for all areas.

If `snakemake` gives 'nothing to be done' after the first iteration, delete the database and PSS/E files and run it again. 
