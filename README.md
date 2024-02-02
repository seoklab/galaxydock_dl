# GalaxyDock-DL (linux only)
GalaxyDock-DL is a protein-ligand docking method which utilizes Conformational Space Annealing(CSA) as a sampling algorithm and deep learning-based scoring functions.

## Installation Guide (linux only)
1. Install [Anaconda](https://www.anaconda.com/products/individual) if you have not installed it yet.<br/>
2. Install [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/download.html) if you have not installed it yet.<br/>
3. Intallation can be done by running below commands in terminal from main directory location. After git clone, below commands should be run in terminal from main directory location.<br/>
4. Clone this Git repository<br/>

```bash
$ git clone git@github.com:seoklab/galaxydock_dl.git
```

Below commands should be run in terminal from main directory location.<br/>

5. Create a conda environment using a env.yml file.<br/>

```bash
$ conda env create --file env.yml
```

6. Activate the conda environment.<br/>

```bash
$ conda activate gd_dl
```

7. Install source files (gd_dl)<br/>

```bash
$ pip install -e .
```

When you use GalaxyDock-DL, please make sure you activate the enviroment in terminal first ("conda activate gd_dl").<br/>

## Usage (linux only)
Below commands should be run in terminal from main directory location.<br/>
We recommend checking src/gd_dl/path_setting.py if you want to change path settings.<br/>

### Running docking
Run docking for a single ligand mol2 file and a protein receptor file without the ligand. (A center coordinate of a docking box (22.5 angstrom^3) is usually set to a coordinate of cognate ligand's center of mass for docking box to include a binding site.)<br/>

Random seeds were set to zero, but you can modify random seeds by adding argument '--random_seed <random_seed value>'.<br/>
Default output directory is set to current working directory, but you can modify random seeds by adding argument '--out_dir <location of output directory>'.<br/>
You can change length of docking box by adding argument '--box_size <box size value in angstrom>'.<br/>

```bash
$ python scripts/run_gd_dl.py -p <path to protein receptor file(.pdb)> -l <path to ligand file(.mol2)> -x <center x coordinate of a docking box> -y <center y coordinate of a docking box> -z <center z coordinate of a docking box>
```

If you want to run docking in terminal from a different directory, you can use bash command with '-d <location of main directory>' below<br/>
```bash
$ python scripts/run_gd_dl_from_other_directory -d <Path to main directory> -p <path to protein receptor file(.pdb)> -l <path to ligand file(.mol2)> -x <center x coordinate of a docking box> -y <center y coordinate of a docking box> -z <center z coordinate of a docking box>
```

### Output files
GalaxyDock_fb.mol2: final output ligand poses sorted by total score<br/>
GalaxyDock_fb.E.info: scores of final output ligand poses in the final bank sorted by total score<br/>
For GalaxyDock_fb.E.info, second column is ranking scores of output poses inferenced by neural network scoring functions. You can ignore values in l_RMSD column since they just represent RMSD calculated by Hungarian algorithm between processed input ligand poses and output ligand poses.<br/>

GalaxyDock_ib.mol2: Initial ligand conformations in the first bank<br/>
box.pdb: Representation of docking box<br/>
GalaxyDock_cl.mol2: clustered final output ligand poses sorted by total score<br/>

Other output files are used during initialization or sampling and not important after docking is finished.<br/>

You can view ligand conformations directly using UCSF chimera<br/>
For example,
```bash
$ chimera GalaxyDock_fb.mol2
```

or you can view ligand conformations and protein receptor
```bash
$ chimera GalaxyDock_fb.mol2 <path to protein receptor file(.pdb)>
```

### Running docking to test

Running docking for the CASF-2016 core set and PoseBusters test set to reproduce the result using SLURM (Simple Linux Utility for Resource Management).<br/>
You can tailor this scripts to work well in your SLURM settings or other linux job schedulers.<br/>
Data directory "total_data/" should be downloaded from [link]() and unzipped into a folder in the main directory.<br/>

This is for CASF-2016 core set
```bash
$ python scripts/multi_run_gd_dl.py test prep
$ python scripts/multi_run_gd_dl.py test run
$ python scripts/multi_run_gd_dl.py test rmsd
$ python scripts/multi_run_gd_dl.py test result
```

This is for CASF-2016 core set using generated molecules generated by CORINA
```bash
$ python scripts/multi_run_gd_dl_corina.py test prep
$ python scripts/multi_run_gd_dl_corina.py test run
$ python scripts/multi_run_gd_dl_corina.py test rmsd
$ python scripts/multi_run_gd_dl_corina.py test result
```

This is for PoseBusters set
```bash
$ python scripts/posebuster_multi_run_gd_dl.py posebuster prep
$ python scripts/posebuster_multi_run_gd_dl.py posebuster run
$ python scripts/posebuster_multi_run_gd_dl.py posebuster rmsd
$ python scripts/posebuster_multi_run_gd_dl.py posebuster result
```

This is for PoseBusters set using generated molecules generated by CORINA
```bash
$ python scripts/posebuster_multi_run_gd_dl.py posebuster_corina prep
$ python scripts/posebuster_multi_run_gd_dl.py posebuster_corina run
$ python scripts/posebuster_multi_run_gd_dl.py posebuster_corina rmsd
$ python scripts/posebuster_multi_run_gd_dl.py posebuster_corina result
```

# License

All code is licensed under MIT license, the weights of the neural networks are licensed under CC BY 4.0.