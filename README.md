# MATPredict: prediction of MAT loci in fungi

<div style="text-align: center;">
  <img src="misc/logo.jpg" alt="Logo" width="500" style="display: block; margin: 0 auto;"/>
</div>

Framework:
The premise for this should be develop model MAT locus representation for taxonomic groups.
Build and use HMMs for specific taxonomic groups
Syntenic searching with something like [https://github.com/gamcil/cblaster](cblaster). 

Description:
This repository is dedicated to the development of a tool that can predict MAT loci in fungal genomes via a machine learning approach. No tool has attempted to predict MAT loci in fungi. This is important because not only are MAT loci important for the evolution and sexual recombination of fungi, but their status is often associated with pathogenicity phenotypes in fungi. Understanding the distribution of MAT loci may facilitate pathogen identification and disease management.

### Preliminary goals: 

Acheive accurate predicton of MAT locus coordinates for every fungal class.

  - This would involve making seperate models for every fungal class, and being able to specify which model to use via command line for users. 

  - Multiple models will be tested to determine which is the most accurate/sensitive. 

Be able to distinguish MAT1-1 from MAT1-2

  - Could I train seperate models using MAT1-1 and MAT1-2?

Integrate prediction accuracy into the logfile

  - Have some kind of % confidence measure 

Be able to show % completion in the logfile

  - 25 percent done...
  - 50 percent done...

Lastly, I want to make this package fast, and make sure the download doesn't take forever. 


### Planned directory structure:
```
MATPredict/
│
├── MATPredict/
│   ├── __init__.py 
│   ├── data_processing.py
│   ├── prediction.py
│   ├── main.py
│   └── utils.py
│
├── tests/
│   ├── __init__.py 
│   ├── test_data_processing.py
│   ├── test_prediction.py
│   ├── test_main.py
│   └── test_utils.py
│
├── model_build/
│   ├── __init__.py
│   ├── get_MAT.py
│   ├── model_build.py
│   ├── get_testset.py
│   └── viz.py
│
├── examples/
│   └── example_prediction.py
│
├── setup.py
├── README.md
├── requirements.txt
└── .gitignore
```
