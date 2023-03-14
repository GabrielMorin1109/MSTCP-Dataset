
# check for change in the dataset
onstart:
    print("Check for database update.")
    shell("Bash/sync_MSWEP.bash")

rule targets:
    input:
        "Data/Data_IBTrACS/ibtracs.ALL.list.v04r00.csv"

rule download_IBTrACS:
    input:
        "Bash/download_IBTrACS.bash"
    output:
        "Data/Data_IBTrACS/ibtracs.ALL.list.v04r00.csv"
    shell:
        """
        {input} {output}
        """

