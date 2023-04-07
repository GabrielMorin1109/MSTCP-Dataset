# TC_Precipitations
A project to practice reproducible analysis while studying Tropical Cyclone Precipitation.

## ENVIRONMENT 
```diff
- TO BE IMPLEMENTED
```
<!-- Create a new environment by running (if you are using mamba) the following command:
```
mamba env create -f environment.yml
``` -->
<!-- ```
conda config --add channels conda-forge
conda config --set channel_priority strict
conda env create -f environment_droplet.yml
``` -->

## DATA 
To access any data, you must accept all provider data licensing. 


- First, to access MSWEP, you have to apply here [http://www.gloh2o.org/mswep/](http://www.gloh2o.org/mswep/). Then, you will receive a Google Drive link containing MSWEP once your request has been approved. You will have to set up `rclone`; you can follow the instructions from gloh2o website here [http://www.gloh2o.org/mswep/#:~:text=What%20is%20the%20most%20effi%C2%ADcient%20way%20to%20down%C2%ADload%20data%20from%20the%20Google%20Drive%3F](http://www.gloh2o.org/mswep/#:~:text=What%20is%20the%20most%20effi%C2%ADcient%20way%20to%20down%C2%ADload%20data%20from%20the%20Google%20Drive%3F). All data will automatically be downloaded from the shared Google Drive. 
- IBTrACS is downloaded automatically from the sever here [https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/](https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/). 

## RUN WORKFLOW
```diff
- TO BE IMPLEMENTED
```

<!-- In order to run the code, download the github, open the file in terminal and activate conda environment with:
```
mamba activate TC_Precipitations
```

Then, run the data analysis workflow with:
```
snakemake --cores 1
```

If there is any change in the yml file, then run :
```
mamba env update -n TC_Precipitations --file environment.yml
``` -->
<!-- 
If necessary, to remove environment, run :
```
mamba remove -n TC_Precipitations --all 
```
-->