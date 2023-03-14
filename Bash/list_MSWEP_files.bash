#!/usr/bin/env bash
# remove old MSWEP_FileName if exist
[ -e Data/MSWEP_FileName.txt ] && rm Data/MSWEP_FileName.txt

# create MSWEP_FileName.txt to echo all filename from GoogleDrive
echo " Size FileName" > Data/MSWEP_FileName.txt

# Loop over NRT and Past folder
for folder in NRT Past
do
    # list all file; select only filename with .nc; with ">>", append filename into MSWEP_FileName.txt 
    rclone ls -v --drive-shared-with-me "GoogleDrive:/MSWEP_V280/$folder/Monthly" | grep ".nc" >> MSWEP_FileName.txt
done