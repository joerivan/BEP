## Directories and documents

* `run_poscar`

	Bash script to submit `build_poscars.py` using slurm. 

* `build_poscars.py`

	Python script to produce POSCAR files based on `data/param.txt` and `data/POSCAR`.

* `run`

	Bash script that starts the VASP and lobster calculations using slurm array jobs.

* `run_vis`

	Bash script to submit `visualise.py` using slurm.

* `visualise.py`

	Python script to produce pDOS and COHP plots based on output files from VASP and lobster calculations.

* data
    + `POSCAR`
	
		POSCAR file with CO on a surface which is geometry optimzed.
	
	+ `param.txt`
	
		Paramter file used to specifiy the max distance to increase the CO molecule in amstrong. This is the first row. The second row is used to specify how many steps are taken.
	
	+ `INCAR`
	+ `KPOINTS`
	+ `POTCAR`
	+ `lobsterin`
* output
	+ images
	
		Directory which will contain the produced images.
	
	+ positions
	
		Directory which will contain the CONTCAR files.

## Dependencies

To make sure everything runs correctly, a local python environment needs to be made. 
In this environment the python packages `ase`, `numpy` and `Pillow` need to be installed.

## Usage

A user should first make sure that the files `data/POSCAR`, `data/param.txt`, `data/INCAR` ,`data/KPOINTS`, `data/POTCAR` and `data/lobsterin` are correct. 
The number of array jobs in `run` should be manually changed.
Then all paths in `run_poscar` and `run_vis` need to be checked.
In these two files the python environment needs to be activated. 
The user should thus make sure that the paths are correctly set in these two files.
When checking the paths in `run_vis`, the user should also check if the range of the for loop correct is. 

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



