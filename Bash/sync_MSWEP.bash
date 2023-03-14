#!/usr/bin/env bash

# create folder Data_MSWEP if doesn't exist
[ -d Data/Data_MSWEP ] || mkdir Data/Data_MSWEP

# syncronize googledrive with local system
for folder in NRT Past
do
    # from GoogleDrive, check if file matching with file in Data_MSWEP
    rclone sync -v --drive-shared-with-me "GoogleDrive:/MSWEP_V280/$folder/Monthly/" Data/Data_MSWEP
done