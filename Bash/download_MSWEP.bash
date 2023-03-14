#!/usr/bin/env bash
# syncronize googledrive with local system
for folder in NRT Past
do
    # from GoogleDrive, check if file matching with file in Data_MSWEP
    rclone sync -v --drive-shared-with-me "GoogleDrive:/MSWEP_V280/$folder/Monthly/" Data/Data_MSWEP
done