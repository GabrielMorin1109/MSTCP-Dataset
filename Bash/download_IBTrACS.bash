#!/usr/bin/env bash
# make the first argument when calling the script as file_path variable
file_path = Data/Data_IBTrACS/

# get IBTrACS and save it to file_dir
wget \
    {$file_dir} \
    https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.ACTIVE.list.v04r00.csv

# create MSWEP_FileName.txt to save 
echo "FileName" > Data/Data_IBTrACS/MSWEP_FileName.txt

# list all file, select only filename with .nc, append filename into MSWEP_FileName.txt 
rclone ls -v --drive-shared-with-me "GoogleDrive:/MSWEP_V280/Past/3hourly" | grep ".nc" > Data/Data_IBTrACS/MSWEP_FileName.txt