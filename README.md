# NMR Data Processing Scripts README

Welcome to the README for a set of scripts designed to facilitate the processing and analysis of Nuclear Magnetic Resonance (NMR) data for molecular modeling. This suite of scripts aids researchers in extracting essential information from the Biological Magnetic Resonance Bank (BMRB) database, converting files into CHARMM-compatible formats, generating restraints, and performing analyses. Below is an overview of each script and its usage.

## `step00_ckbmrb.sh` (BMRB Data Retrieval Script)

### Overview

This script retrieves diverse NMR-related data from the BMRB database for a given model (PDBID). It validates the presence of general distance restraints, hydrogen bond data, dihedral angles, and dipolar coupling information.

### Usage

1. Ensure that you have the `csh` (C Shell) interpreter installed.
2. Run the script with the model identifier (PDBID):


Example: `./step00_ckbmrb.sh 1A7F`

3. The script fetches and assesses the availability of various data types.
4. Results are appended to the `NMR_ID.sort.ckbmrb.txt` file.

## `step01_pdb2cpdb.sh` (PDB to CHARMM PDB Format Conversion Script)

### Overview

This script converts Protein Data Bank (PDB) formatted files into the CHARMM PDB format, which is crucial for molecular dynamics simulations and analysis using the CHARMM software suite.

### Usage


Example: `./step01_pdb2cpdb.sh 1A7F.pdb`

## `step02_getrestraint.sh` (NMR Restraint Data Downloader and Processor)

### Overview

This script simplifies the download and processing of NMR restraint data from the BMRB server. It streamlines the preparation of NMR data for further analysis and modeling. The script automates the retrieval of different restraint types, including general distance restraints and dihedral angle restraints.

### Usage


Or, using a previously generated file: NMR_ID.sort.ckbmrb.txt


Example: `./step02_getrestraint.sh 1A7F 1 0 0 0`

## `step03_mkrsr.sh` (Create CHARMM Readable Restraints)

### Overview

This script assists in preparing and generating restraints (e.g., NOE, dihedral) for CHARMM simulations. It manages restraint-related files and directories, while also checking flags and existing files to track the restraint generation process. Please note that additional context is necessary to fully comprehend the script due to commented sections and potential variations in usage.

### Usage


Example: `./step03_mkrsr.sh 1A7F.pdb`

## `step04_ckaqua.sh`

### Overview

This script conducts checks and analyses relevant to CHARMM and AQUA restraints in molecular modeling. The script involves tasks like extracting file names, configuring directory paths, and executing processes based on specific conditions.

### Usage


Example: `./step04_ckaqua.sh 1A7F.pdb`

## `step05_runref.sh` (Generate CHARMM TrioSA Refinement and Calculate Restraint Violations)

### Overview

This script facilitates refining models using the CHARMM TrioSA method and subsequently calculating restraint violations (e.g., NOE, dihedral) for each model in the ensemble.

### Usage


Example: `./step05_runref.sh 1A7F.pdb`

# Issue Uploading Source Files for Review Submission


It appears that the size limitations set by the submission platform have prevented me from successfully including the convpl, STAP, and CHARMM files alongside my manuscript. Unfortunately, these files are integral to a comprehensive understanding of my research.

I'm writing to kindly request guidance on the best course of action in this situation. It's of utmost importance to me that I adhere to the journal's guidelines and ensure the successful review of my work.

If there's an alternative method to submit these files, or if I need to modify them in any way, please do let me know. I'm committed to working together to find a solution that aligns with the journal's requirements and facilitates a thorough review process.

For any further information or clarification, feel free to reach out to me directly at yycho@kribb.re.kr. Your prompt attention to this matter would be greatly appreciated.

Feel free to use and adapt these scripts as needed for your research purposes. If you have any questions or need further assistance, please don't hesitate to reach out.

