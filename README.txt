**************************************************************************

THESIS DATA DIRECTORY 


This directory contains test files, the dataset and results generated 
during the 2021/22 SAAPdap/pred MSci project reruns done by Barbara Mikucka. 

**************************************************************************

1. File: Q6GZV6.txt
	SwissProt file used to test sprotFTdist.py preformance during 
	development (when input used was SwissProt file instead of PDB file)

2. Dir: HumVar_data
	The directory contains 2 files used for SAAPdap runs. 
	PD and SNP HUmVar 2011/12 datasets with protein UniProt acccesion code 
	and mutation info.

	2. 1. humvar-2011_12.deleterious.pph.txt
	2. 2. humvar-2011_12.neutral.pph.txt

3. Dir: JSON_SNP
	The directory contains all SNP JSON files - outputs of SAAPdap runs 
	using the SNP HumVar dataset. (21,151 files)
	Details on the runs: https://github.com/bmikucka/project

4. Dir: JSON_PD
	The directory contains all PD JSON files - outputs of SAAPdap runs 
	using the PD HumVar dataset. (22,196 files)
	Details on the runs: https://github.com/bmikucka/project

5. File: snp.csv
	Includes all SNP datapoints used for SAAPpred. Contains information 
	from all non-error SNP JSON files from SAAPdap.

6. File: pd.csv
	Includes all PD datapoints used for SAAPpred. Contains information from 
	all non-error PD JSON files from SAAPdap.

7. File: inputs.dat  
	Contains list of features to include in Weka Machine Learning analysis 
	- excludes SprotFTdist results. 

8. File: inputs_updated.dat
	Contains list of features to include in Weka Machine Learning analysis 
	- includes SprotFTdist results.

9. Dir: PRED_MCC
	Contains all MCC result files for all runs conducted. Results included 
	were taken from Weka output files. 
	File structure: a line for each Weka ML run: 
	(1) test run number (2) testing MCC value (NOT training MCC).
	Run 1 is not included - script errors prevented generating results.
	Details on run set up is provided in MCC_results.xlsx
	More details: https://github.com/bmikucka/project

10. File: MCC_results.xlsx
	Details on parameters used for each Weka Machine Learning run.

**************************************************************************
