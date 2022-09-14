# What this app does:
This is a web app containing exploratory graphs for spatialomics data. 
It has been designed while working with
spatially resolved proteomics and transcriptomic data from the GeoMx DSP 
and will draw
5 interactive plots on any dataset that you give it. The plots are:
* A parallel categories plot of the selected columns
* A scatter plot of the data (with spearman's rho and p-value)
* 3D UMAP
* Volcano plot of the data (with univariate results of a Linear Mixed Model)
* And a univariate Cox regression horizontal bar plot

# Notes on data input
While the app is robust in reading file formats, it assumes that the uploaded data
 is formatted so that all information is stored in columns. All biomarker counts
 need to be in columns and all factors also need to be in columns and one column
 needs to be some kind of identification to track the source (patient). 

# Quick start:
The main script is app/mk_dash.py but before running it you need to setup an environment,
this can be done using [Conda](https://conda.io/docs/user-guide/install/index.html).

## Windows users:
* After installing find the anaconda prompt
* Use anaconda prompt to run the specified commands!

### Linux and mac users:
* You probably know how to use conda and a terminal! :) 

_(if you don't: look for the link to conda docs above!)_

## Open a terminal in path/to/mk_dsp_explorer and run:
```conda env create -n dashenv -f env/environment.yml```

## Next, activate your new environment:
```conda activate dashenv```

* In the future you may wish to list your environments:
```conda env list```

## Configuration 
All settings made in the app are now also written to conf.json, this file can be
backed up when switching. 
## When you are ready you can run the main app by doing:
```python app/mk_dash.py```

# NOTES:

Expected data types of the different analyses:
* Parallel plot - categorical data will be transformed to numbers ordered by sorting
otherwise colors won't work!
* LMM - binary or ordinal encoding of the fixed effect, continuous values for biomarkers
* Cox - binary, ordinal or continuous - estimates a ratio that is influenced by the order of values
* Correlation - continuous (the idea is to check for correlation of various biomarkers)
* 3D UMAP - distance is measured for the biomarkers (continuous logscale) could be improved
to not be so sensitive about strings for coloring clusters.

# Acknowledgements:
Some really top notch software that go into making this app work are:

[Conda](https://conda.io/)

[Lifelines for Python](https://lifelines.readthedocs.io/en/latest/)

[Dash](https://dash.plot.ly/dash-core-components)

[Dash Bootstrap Components](https://dash-bootstrap-components.readthedocs.io/en/latest/)

[Plotly Dash](https://plotly.com/dash/)

[Plotly](https://plotly.com/)

[Numpy](https://numpy.org/)

[Pandas](https://pandas.pydata.org/)

[Statsmodels](https://www.statsmodels.org/)

[scikit-learn](https://scikit-learn.org/)



