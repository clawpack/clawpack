### Folder Organization
* **adjoint:**

Contains code to solve the adjoint problem.

The output times specified in this directory should agree with those for the forward code.

### Running the Code

* Go to the folder **adjoint** and run in a terminal:

```
python maketopo.py
make new
make .plots
```

The code will produce two new folders: _output and _plots. 
The first one contains all the output files, while the latter one contains the plots and interactive 
visualization apps.

* Go to the main folder **tsunami_Alaska** and run in the terminal:

```
python maketopo.py
make new
make .plots
```

The code will produce two new folders: _output and _plots. 
The first one contains all the output files, while the latter one contains the plots and interactive 
visualization apps. Note that this simulation takes a long time to complete.

### Running Variations

* Running the example with adjoint flagging:

Run in the terminal:

```
python run_adjoint_flagging.py
```

This example can be run with either the regular surface-flagging technique (by using the flag2refine file) or with error-flagging (by using errf1 file). These options can be set by setting the flag2refine and/or flag_richardson flags in setrun.py.