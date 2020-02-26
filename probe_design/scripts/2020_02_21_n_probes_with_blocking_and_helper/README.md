# smFISH Probe Design
Probe design pipeline for plasmid or transcript imaging experiments.

## Overview
This pipeline enables design of complex oligo probe sets targeted at mobile genetic elements or transcripts. The main pipeline is a snakemake workflow.

## Before running the pipeline
1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html),
2. Install the environment by either
   - [IGNORE THIS] Running `conda env create -f smfish.yml` in a Terminal window,\

  OR

   - Running the following command\
     `conda env create smfish python=3.5`\
     `conda install pandas=0.25`\
     `conda install -c anaconda biopython`\
     `conda install -c etetoolkit ete3`\
     `conda install snakemake`\
     `conda install blast`\
     `source activate smfish`\
     `pip install SetCoverPy`
     `pip install tables`
     `pip install openpyxl`
     `pip install matplotlib`
     

3. Activate the environment by running `source activate smfish`,
4. Edit the `smfish_config.json file` to point the pipeline to the correct directories.
   - `__default__`
      * `SCRIPTS_PATH`: path to the folder that contains all the scripts
      * `DATA_DIR`: path to the folder that contains input folders and files
   - `blast`
      * `custom_blast`: path to the custom reference fasta file. All probes will be blasted against this fasta file, in addition to features that are contained in the plasmid or transcripts themselves.
   - `params`
      * `length_threshold`: length threshold for features to be analyzed. Any feature shorter than this length will be ignored in the design. Default is 50.
      * `theme_color`: color scheme for probe design summary plots. To use figures in a white background, use "black"; for figures in a black background, use "white" or any other colors that are not black. 
   - `simulations`
      * `simulation_table`: path to the simulation summary file

## Input
1. Simulation summary file (simulation_table_test.csv)
   - A csv file containing all the designs to be run.
      * `DESIGN_ID`: identifier for each design
      * `SAMPLE`: name of the input FASTA file without file extension
      * `MAX_CONTINUOUS_HOMOLOGY`: maximum continuous homology (measured in bp) for a probe-target hit to be considered significant. Lower values leads to more stringent designs. Default is 14 bp.
      * `MIN_TM`: minimum melting temperature threhold
      * `MAX_TM`: maximum melting temperature threhold
      * `GC`: minimum probe GC content threhold
      * `INCLUDE_START`: number of nucleotides to exclude at the beginning of the feature sequence
      * `INCLUDE_END`: number of nucleotides to exclude at the end of the feature sequences
      * `PROBE_SELECTION_METHOD`: method for selecting probes. Available options are
         1. `SingleBestProbe`: select the top probe for each feature, if available
         2. `AllSpecific`: select all probes that are specific and only specific to its target feature
         3. [NOT IMPLEMENTED YET]`AllSpecificPStartGroup`: select all probes that are specific and only specific to its target feature within each segment of the feature sequences. By default the feature sequences are dividied into block resolutions of 100bp regions. If there are less than 15 probes available (average one probe per block), the block resolution is modified in 20bp decrements until there are 15 probes or the block resolution is zero, whichever happens first.
         4. `MinOverlap`: select all probes that are specific and only specific to its target feature with minimum overlap in their target coverage
         5. `TopN`: select the top *n* probes for each feature
       * `PRIMERSET`: primer sets to include in the final probes. There are three sets (A, B, and C) availble in the current version. User specific primer sets can also be added if necessary.
       * `TPN`: number of top probes to select for each taxon, if the probe selection method is set to `TopN`
       * `BOT`: minimum blast on target rate threshold. Probes with blast on target values lower than this value is considered *promiscuous*, and is not included in the final probe pool.
       * `BARCODESELECTION`: method for barcode assignment to taxa. Available options are:
         1. MostSimple: assign barcodes by barcode complexity, starting with the simplest ones. Barcodes with more bits are considered more complex.
         2. Random: randomly assign barcodes to taxa
         3. MostComplex: assign barcodes by barcode complexity, starting with the most complex ones. Barcodes with more bits are considered more complex.
       * `BPLC`: minimum blocking probe length threhold. Blocking probes with length lower than this threshold is considered likely to be washed off and do not need to be included in the final probe pool. Default is 15 bp.
       * `BT`: maximum bit score threhold for off target blast hits. Probes with maximum off target hit bit scores larger than this threshold will be removed from the complex pool. Default is 25.
2. Plasmid file
   - A snapgene .dna file containing plasmid sequences with annotations, placed at [DATA_DIR]/[SAMPLE]/input.
3. Custom blast database file
   - A FASTA file containing the reference sequences. Typical examples include whole genome sequences of the organisms that will be used in the experiments, also placed at [DATA_DIR]/[SAMPLE]/input.

## Output

1. Simulation results file
   - A csv file containing all the parameters for all the designs, as well as some summary statistics for each design
2. Probe folder
   - A folder containing selected probe summary files for each taxa, a concatenated file containing all selected probes, a file containing information for all the blocking probes, as well as text files that can be sent as is to array synthesis vendors for complex oligo pool synthesis.

## Running the pipeline
Run `snakemake --configfile smfish_config.json -j n`, where `n` is the number of cores to be used. If the pipeline excuted without errors, you should see a file called `simulation_table_test_results.csv` in the same directory where you put the `simulation_table_test.csv` file. It can be useful to run a design at a high taxonomic rank (phylum, for example) to make sure that the pipeline runs correctly with the input files.
