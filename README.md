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
Run docking for a single ligand mol2 file and a protein receptor file without the ligand. (A center coordinate of a docking box (22.5 angstrom^3) is usually set to a coordinate of cognate ligand's geometric center for docking box to include a binding site.)<br/>

Random seeds were set to zero, but you can modify random seeds by adding argument
```bash
--random_seed <random_seed value>
```

Default output directory is set to current working directory, but you can modify random seeds by adding argument
```bash
--out_dir <location of output directory>
```

You can change length of docking box by adding argument
```bash
--box_size <box size value in angstrom>
```

```bash
$ python scripts/run_gd_dl.py -p <path to protein receptor file(.pdb)> -l <path to ligand file(.mol2)> -x <center x coordinate of a docking box> -y <center y coordinate of a docking box> -z <center z coordinate of a docking box>
```

Example for 3rsx

```bash
$ python scripts/run_gd_dl.py -p ./example/3rsx/3rsx_protein.pdb -l example/3rsx/3rsx_ligand.mol2 -x 69.637 -y 49.989 -z 10.160 --out_dir example/output_dir/
```

If you want to run docking in terminal from a different directory, you can use bash command with '-d <location of main directory>' below<br/>
```bash
$ python scripts/run_gd_dl_from_other_directory -d <Path to main directory> -p <path to protein receptor file(.pdb)> -l <path to ligand file(.mol2)> -x <center x coordinate of a docking box> -y <center y coordinate of a docking box> -z <center z coordinate of a docking box>
```

### Output files
- `GalaxyDock_fb.mol2`: Contains the final output ligand poses, sorted by total score.
- `GalaxyDock_fb.E.info`: Provides the scores of the final output ligand poses in the final bank, sorted by total score.

For `GalaxyDock_fb.E.info`:
- The second column (Energy) shows the ranking scores of output poses inferred by neural network scoring functions.
- You can ignore the values in the l_RMSD column, as they only represent RMSD calculated by the Hungarian algorithm between processed input ligand poses and output ligand poses.
- You can also ignore the other columns, which correspond to the values of GalaxyDock BP2 Score energy components multiplied by their weights (ATDK_E: AutoDock Energy, INT_E: AutoDock intra-ligand energy, DS_E: Drug Score, HM_E: Hydrophobic interaction, PLP: PLP score).

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
Data directory "total_data/" should be downloaded from [link](https://drive.google.com/file/d/1sILlVoda3_f6E3oc4Rr8Jn6Ae93zA7q4/view?usp=sharing) and unzipped into a folder in the main directory.<br/>

To run the scripts without SLURM, prefix the script filename with `noslurm_`. Since docking tests can be time-consuming, you can take advantage of multi-CPU support in **run mode** provided by the `noslurm_*` scripts. For instance, you can execute the following command to utilize multiple CPUs:

```bash
python scripts/noslurm_multi_run_gd_dl.py run {n_cpu}
```

This is for CASF-2016 core set
```bash
$ python scripts/multi_run_gd_dl.py prep
$ python scripts/multi_run_gd_dl.py run
$ python scripts/multi_run_gd_dl.py rmsd
$ python scripts/multi_run_gd_dl.py result
```

This is for CASF-2016 core set using generated molecules generated by CORINA
```bash
$ python scripts/multi_run_gd_dl_corina.py prep
$ python scripts/multi_run_gd_dl_corina.py run
$ python scripts/multi_run_gd_dl_corina.py rmsd
$ python scripts/multi_run_gd_dl_corina.py result
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

### Precompiled Binary

We recommend using the precompiled 'ligdock' binary file located in `src/gd_dl/bin/`.

### Compiling from Source

If you wish to compile it yourself, navigate to the `binary_src/` directory by running `cd binary_src/` and then execute the following commands (note the ".." at the end of the second command):

```bash
$ mkdir -p build && cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort ..
$ cmake --build . --target all -j8
```

This will generate a new 'ligdock' binary file in the `binary_src/bin/` directory.

Ensure that the following compilers are installed before compiling:

- **icc (Intel C Compiler)**

```
icc (ICC) 2021.9.0 20230302
Copyright (C) 1985-2023 Intel Corporation.  All rights reserved.
```

- **ifort (Intel Fortran Compiler)**

```
ifort (IFORT) 2021.9.0 20230302
Copyright (C) 1985-2023 Intel Corporation.  All rights reserved.
```

If you want to replace the original binary file with the new one, you can execute the following command:

```bash
$ cp ../bin/ligdock ../../src/gd_dl/bin/ligdock
```

# Citation

If you utilize this code or the models in your research, please cite the following paper:
```
@article{lee2024galaxydock,
  title={GalaxyDock-DL: Protein--Ligand Docking by Global Optimization and Neural Network Energy},
  author={Lee, Changsoo and Won, Jonghun and Ryu, Seongok and Yang, Jinsol and Jung, Nuri and Park, Hahnbeom and Seok, Chaok},
  journal={Journal of Chemical Theory and Computation},
  year={2024},
  publisher={ACS Publications}
}
```

# License

All code, except for the code in the "binary_src" directory, is licensed under the MIT license. The weights of the neural networks are licensed under the CC BY-NC 4.0 license. Code in the "binary_src" directory is licensed under the CC BY-NC-ND 4.0 license.
