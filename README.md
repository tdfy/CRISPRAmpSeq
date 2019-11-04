# CRISPRAmpSeq
CRISPR edit analysis by amplicon sequencing 

```
Configuration Table Example:
 Directory                       File_Name               Target      Sample  Design    Experiment_Name      Time_Point    Chron_Order  Noise
------------------------------  ----------------------  --------  --------  --------  -----------------  ------------  -------------  -------
/mnt/c/Export/Amp_test/PACE_21  90930_B2M_1.output.csv  B2M              1  T         PACE_20_B                     1              1  Y
/mnt/c/Export/Amp_test/PACE_21  90930_B2M_2.output.csv  B2M              2  C         PACE_20_B                     2              2  Y

usage: CRISPRAmpSeq.py [-h] -p PATH -f FILE

Application to call CRISPR edits from amplicon sequencing libraries

optional arguments:
  -h, --help
  -p PATH, --path PATH  path to configuration file...
  -f FILE, --file FILE  name of configuration file.[CSV]
```

