# NoiseXcor

This is a MATLAB code that will compute noise cross correlations. 

---
### UPDATE 

As of 1 July 2021 the __master__ branch has been renamed __main__. You need to execute the following code to make this change in your local repository.

```  git branch -m master main
git fetch origin
git branch -u origin/main main
git remote set-head origin -a
```

### Requirements

1) The _NoiseXcor_ package relies on the [GISMOTOOLS](https://geoscience-community-codes.github.io/GISMO/) _waveform_ object to hold the data and metadata for each  correlation waveform. You can follow the installation instructions [here](https://github.com/geoscience-community-codes/GISMO/wiki/Getting-Started). You will want to follow all steps so that the GISMOTOOLS library is loaded every time you start MATALB.

The remainder of the code uses internal MATLAB libraries. The current version of the code depends on MATLAB v9.2 (2017a). We are constantly working to utilize internal improvements in MATLAB.

### Installation

After cloning this repository to your local machine, you will need to add the _src_ folder to your MATLAB path. Below is an example in which we append the _src_ folder to the end of the MATLAB path. 

```
addpath('/Users/dmikesell/GIT/NoiseXcor/src', '-end');
```

You may want to permanently add this path so you do not have to run this command each time you start MATLAB. You can do this by running the path GUI.

```
pathtool
```

### Example

See the tutorial in _examples/Example\_SAC/README.md_ to see how the code runs.
