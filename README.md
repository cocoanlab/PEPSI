# A Computational Mechanism of Cue-Stimulus Integration for Pain in the Brain

[https://doi.org/10.5061/dryad.41ns1rnpj](https://doi.org/10.5061/dryad.41ns1rnpj)

The codes and data for the manuscipt,

"A Computational Mechanism of Cue-Stimulus Integration for Pain in the Brain"

PEPSI stands for **P**ain **E**xpectation and **P**ain **S**timulus **I**ntegration.

## Description of the data and file structure

### codes

##### step1\_CalcSubspacesEncodingperf.m

* Implements the processes and visualizations depicted in Figs. 2A-B (Calculation of subspaces and encoding performances).

#### step2\_VisTraj.m

* Implements the visualization of the trajectories in subspaces.

##### step3\_Traj2Behv.m

* Implements the processes and visualizations depicted in Fig. 5 (Reconstructing behavioral patterns from the neural trajectories).

### data

##### data\_for\_replication.mat

* Contains all data necessary for implementing the steps in the codes folder.
  * `neurAvg.FIR`       : subject averaged FIR response. comprising whole voxel.
  * `neurAvg.CuetBeta`, `neurAvg.StimtBeta`       : temporal encoding weights using `neurAvg.FIR` as Y, and `CueStimX` as X
  * `CueStimX`      : [Intercept, CueInfo, StimInfo]. Normalized. Condition ordered as in      variable "condName"
  * `parcelIndx`: spatial index of each network
  * `behvout`: mat size of #cond X #sub. pain reports. condition is ordered as in variable "condName"
  * `cmaps`: color maps for "condName"
  * `cmaps7`: color maps for 7 large-scale functional networks. 

##### template.nii

* Provides the brain template for variable "neurAvg" in the data_for_replication.mat.
  * Example usage
    * `obj = fmri_data(fullfile(basedir, 'data', 'template.nii'))`. This will results in same voxel size in data of the `neurAvg`.

## Sharing/Access information

Links to other publicly accessible locations of the data:

* [https://github.com/cocoanlab/PEPSI](http://...)

## Code/Software

Codes tested on Ubuntu 20.04, matlab R2021b.

Dependencies are...

* SPM12
* [https://github.com/cocoanlab/cocoanCORE.git](https://github.com/cocoanlab/cocoanCORE.git)
* [https://github.com/canlab/CanlabCore.git](https://github.com/canlab/CanlabCore.git)
* [https://github.com/didch1789/yanchogosu_toolbox.git](https://github.com/didch1789/yanchogosu_toolbox.git)

