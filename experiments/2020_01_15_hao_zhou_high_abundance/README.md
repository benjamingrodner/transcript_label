# README for 2020_01_15 high abundance transcript in mouse gut experiment
Based on RNA seq data paired with metagenomic sequencing done by Hao Zhou in Brito lab.
Sequence data in ../../probe_design/data
Probe design in ../../probe_design/runs/DSGN_001
Probes ordered in ../../orders/2020_01_07_high_abundance_transcripts.xlsx

## Targeted Transcripts
- ndhM:
    - aerobic-type carbon monoxide dehydrogenase, large subunit CoxL
    - d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__UBA9475;s__
    - rpkm_average = 4,648
- argC:
    - Belongs to the NAGSA dehydrogenase family. Type 1 subfamily
    - d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__DTU089;g__Acutalibacter;s__Acutalibacter muris
    - rpkm_average = 3284
## Probes
Transcript probes:
p.ndhM_322.27.27
p.ndhM_302.27.27
p.argC_553.25.25
p.argC_529.25.25

rRNA 16s probe:
EUB_RRX

Branch Probes:
b.27.*27.*27.*28
b.25.*25.*25.*30

Adapter probes:
a.28.*28.*R1
a.30.*30.*R7

Readout probes:
R7* (Alexa 647)
R1* (Alexa 488)

## Results
EUB stain failed. Some spots from transcripts, but I could not find bacteria in the sections
