#!/usr/bin/env bash

# file_path variable defining the path where IBTrACS data will be stored
file_path="Data/Data_IBTrACS/"

# create folder Data_IBTrACS if doesn't exist
[ -d $file_path ] || mkdir $file_path

# get IBTrACS and save it to file_path
wget -O ${file_path}/ibtracs.ALL.list.v04r00.csv \
    https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.ALL.list.v04r00.csv