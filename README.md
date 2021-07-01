# NoiseXcor

This is a MATLAB code that will compute noise cross correlations. 

---
## UPDATE: 

As of 1 July 2021 the __master__ branch has been renamed __main__. You need to execute the following code to make this change in your local repository.

```  git branch -m master main
git fetch origin
git branch -u origin/main main
git remote set-head origin -a
```

### Preliminaries:

1) The _NoiseXcor_ package relies on the [GISMOTOOLS](https://geoscience-community-codes.github.io/GISMO/) _waveform_ object to hold all data and metadata for each waveform. You can follow the installation instructions [here](https://github.com/geoscience-community-codes/GISMO/wiki/Getting-Started). You will want to follow all steps so that the GISMOTOOLS library is loaded every time you start MATALB.

The remainder of the code uses internal MATLAB libraries. The current version of the code depends on MATLAB v9.2 (2017a). We are constantly working to utilize internal improvements in MATLAB.

See the Example_SAC/README.md and the example to see how the code runs.

#### Authors: Dylan Mikesell and Piero Poli
Send comments to: _dylanmikesell at boisestate dot edu_ and _ppoli at mit dot edu_.
