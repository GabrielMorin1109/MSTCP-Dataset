#!/usr/bin/env bash
# make the first argument when calling the script as folder_path variable
folder_path=$1

# create Data_MSWEP if doesn't exist in specified fodler
[ -d "Data_MSWEP" ] && rm $folder_path/Data_MSWEP

# syncronize googledrive with local system
# for folder in NRT Past
# do
#     # from GoogleDrive, check if file matching with file in Data_MSWEP
#     rclone sync -v --drive-shared-with-me source:"GoogleDrive:/MSWEP_V280/$folder/Monthly" $folder_path/Data_MSWEP
# done

for file_name 202301.nc 202212.nc
do
    rclone sync -v --drive-shared-with-me source:"GoogleDrive:/MSWEP_V280/NRT/Monthly" $folder_path/Data_MSWEP
done
