#!/usr/bin/env bash
# from GoogleDrive ($1 folder), check if file matching with file in Data_MSWEP and update 
sync_from_GoogleDrive_folder() {
    rclone sync -v --drive-shared-with-me "GoogleDrive:/MSWEP_V280/$1/Monthly/" Data/Data_MSWEP/$1
}

# Get the time of last data modification (only keep year, month and day). If files do not exist, it is not a problem
date_last_modif=$(date -r Data/Data_MSWEP/ +'%Y-%m-%d')
# Get current date minus 7 days. It will verify if last update occurred in the past week.
date_of_last_week=$(date -d 'now - 7 days' +'%Y-%m-%d')


# If [ last modification is an older than (-ot) a week old ] OR [ if Data_MSWEP doesn't exist (-e) ], then run code.
if [ $date_last_modif -ot $date_of_last_week ] || [ ! -e "Data/Data_MSWEP/" ]
then
    echo "The last synchronization was over a week ago. Synchronized database locally from GoogleDrive."

    # Recursively create folder Data_MSWEP if doesn't exist (will not overwrite)
    mkdir -p Data/Data_MSWEP/NRT
    mkdir -p Data/Data_MSWEP/Past

    # from GoogleDrive (NRT and Past folder), check if the files are matching with the files in Data_MSWEP
    sync_from_GoogleDrive_folder NRT
    sync_from_GoogleDrive_folder Past
fi

################################
# Bash condition cheat sheet (https://unix.stackexchange.com/a/304258)
#     "A ; B" Run A and then B, regardless of success of A
#     "A && B" Run B if A succeeded
#     "A || B" Run B if A failed
#     "A &" Run A in background.

# Advance Bash-Scripting Guide: 8.1. Operators (https://tldp.org/LDP/abs/html/ops.html#LOGOPS1)