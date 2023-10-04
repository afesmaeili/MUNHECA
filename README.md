## MUNHECA-v1.0
MUNHECA is a public astrophysical code written in python to simulate the electromagnetic cascade and compute the neutrino spectrum from the cascade development. 

If you used MUNHECA please cite arxiv: 2310.01510

## INSTALLATION AND EXECUTION
MUNHECA uses the numpy and Scipy packages and does not need any installation. 

To run the code, a .txt file should be created which contains the input to the code. 
An example of the input file, the test.txt, is provided in the /Work directory.   

To run the Monte Carlo session, the following command can be used 

```
python [PATH/TO/MUNHECA]/run/Main.py /[PATH/TO/INPUT]/[InputFileName].txt
```

After the execution, the results of this session will be saved, together with a copy 
of the input file, to the /results directory. The copied input file's name, depending on 
the chosen injection spectrum (Monochrome) or (PowerLaw) will be 0M_input.txt or 0S_input.txt, respectively.

To obtain the neutrino fluxes at the Earth the relevant syntaxes are:

```
python [PATH/TO/MUNHECA]/run/nuSpec.py ../results/[DESTINATION_DIR]/0M_input.txt
```

or

```
python [PATH/TO/MUNHECA]/run/nuSpec_Weight.py ../results/[DESTINATION_DIR]/0S_input.txt 
```

for (Monochrome) or (PowerLaw) injections, respectively. 

These commands end in creating a text file NEUTRINO_EARTH.txt in the same destination directory, 
which contains a table of the neutrinos fluxes, containing all-flavors and each flavor separately, at the Earth. 


## AUTHORS
- AmirFarzan Esmaeili
- Arman Esmaili
- Pasquale Dario Serpico




