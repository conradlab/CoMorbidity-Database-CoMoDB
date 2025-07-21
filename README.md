# Mapping whole-organism genetic comorbidities across model organisms using unified ontologies

## CoMoDB
This repository contains code and files for our article in preparation in GENETICS (Peaslee et al, 2025).

CoMo DB is an R Studio analysis pipeline that uses an input gene list to then collect and tally phenotypes per biological system in five organisms (human, mouse, zebrafish, fruit fly, and roundworm) from a ![MARRVEL API](marrvel.org).

## Analysis Instructions

### Setup Input Files
CoMo DB requires 1) gene list in txt file format (one gene per line) 2) phenotype ontology terms for each species (column 1 contains phenotype ID's and column 2 contains phenotype term names) 3) R code "MARRVEL_COMO_DB.R"


```
Required files: 
├── mouse_mgi.txt                     # List of mouse NOA genes
├── fertility_categories/            # Folder containing phenotype ontology terms used in "reproduction" MCat
│   ├── human_terms.txt              
│   ├── mouse_terms.txt
│   ├── zebrafish_terms.txt
│   ├── fly_terms.txt
│   └── worm_terms.txt
├── MARRVEL_COMO_DB.R                # R script(s) for preprocessing and analysis

# Optional reference file:
├── fertility_cats_big_table_formatted_june2024.txt
│                                     # Big table of reproduction phenotype ontology terms organized by species.
```

### Run Pipeline

MARRVEL_COMO_DB.R should be opened with R Studio and input files should be loaded into environment with code file.

Output file "big.out" is used for downstream analysis of phenotypes per gene, per species. 

## Contact us
![Conrad Lab Website](conradlab.org)

Caitlin Peaslee (peasleec@ohsu.edu)
Don Conrad (conradon@ohsu.edu)

## Funding Acknowledgement
This work was supported by grants R01HD078641 and P50HD096723 from the National Institutes of Health. D.F.C. is supported in part by the National Institute of Health Office of Directors (NIH/OD) Grant P51 OD011092 to the Oregon National Primate Research Center.
