## Directories and documents

* `run_poscar`
bash script for submitting `build_poscars.py` using slurm. 
* `build_poscars.py`
* `run`
* `run_vis`
* `visualise.py`
* data
    + `POSCAR`
	poscar file
	+ `param.txt`
	+ `INCAR`
	+ `KPOINTS`
	+ `POTCAR`
	+ `lobsterin`
* output
	+ images
	+ positions

## Dependencies

To make sure everything runs correctly, a local python environment needs to be made. 
In this environment the python packages `ase`, `numpy` and `Pillow` need to be installed.

## Usage

A user should first make sure that the files `data/POSCAR`, `data/param.txt`, `data/INCAR` ,`data/KPOINTS`, `data/POTCAR` and `data/lobsterin` are correct. 
Then all paths in `run_poscar` and `run_vis` need to be checked.
In these two files the python environment needs to be activated. 
The user should thus make sure that the paths are correctly set in these two files.

Then the user can start by submitting the `run_poscar` file with slurm.

```bash
sbatch run_poscar
```

After this is done the user can then submit `run`.

```bash
sbatch run
```

Lastly the user should submit the `run_vis` file when the VASP/lobster calculations are done.

```bash
sbatch run_vis
```

Now all calculations should be complete and the images should be inside the `output/images` directory.



