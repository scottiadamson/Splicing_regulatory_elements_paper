# Vex-seq 2.0 Experimental Protocol

A breif schematic describing the plasmid assembly is availble in Figure 2A [here](https://www.biorxiv.org/content/10.1101/2021.05.14.444228v1)

## Oligonucleotide Sequence Design

All sequences used in Vex-seq have the following general layout:  
5’ CTGACTCTCTCTGCCTC test_sequence CAATTGACTACTAGT 10_nt_barcode TCTAGAGGGCCCGTTTA 3’  
Neither the test sequence or barcode should contain MfeI (CAATTG) or SpeI (ACTAGT) restriction enzyme sequences.  

## Plasmid Assembly

1. Digest at least 3 μg of [modified pcAT7-Glo1](https://www.addgene.org/160996/) with PstI and XbaI. Incubate at 37°C for one hour and heat inactivate at 65°C for 20 minutes. Purify digested DNA using column-based purification (such as Zymo DNA Clean & Concentrator).
2. Blunt ends of digested plasmid (product from step 1) with (1μL/μg DNA) DNA Polymerase I, Large (Klenow) Fragment in Cutsmart buffer supplemented with 33 μM of each dNTP. Incubate at 25°C for 15 minutes. Add (1μL/μg DNA) Calf Intestinal Alkaline Phosphatase (CIP) to reaction and incubate at 37°C for 30 minutes. Inactivate reaction with final concentration of 10 mM EDTA and heat inactivate at 75°C for 20 minutes.
3. Electrophorese blunted and digested DNA fragments (product from step 2) on a 0.5% agarose gel. Gel purify large band (~5.6 kb). If purified DNA has high contamination (260/230 < 2.0) after gel purification, use column purification or ethanol precipitation to further remove contaminants.
4. PCR amplify diluted (20 ng/μL) oligonucleotide pool with Oligo_FWD and Oligo_REV for 10 cycles with 52°C annealing temperature. Purify PCR product with column-based purification.
5. Gibson assemble amplified oligonucleotide pool (product of step 4) into blunted and digested vector (product of step 3). Assemble solution containing a 2:1 molar ratio of insert to vector and NEBuilder® HiFi DNA Assembly Master Mix. Use 100 ng of digested vector per reaction and make two or three total reactions. Include an insert negative control in parallel (for steps 5-7) to estimate number of colonies from re-circularized vector. Incubate reaction at 50°C for 1 hour. Purify assembled reaction using ethanol precipitation or column-based purification. Resuspend or elute purified DNA in 10 μL. High purity DNA is required for efficient transformation. 
6. Electroporate assembled DNA (from step 5) 2.5 μL at a time in electrocompetent E. coli (i.e. E. cloni® 10G Electrocompetent Cells from Lucigen). Follow manufacturer’s protocol for recovery in media. Pool media from all transformations after growth. Use serial dilutions to make a 1:10,000 dilution and plate 100 μL on standard sized LB ampicillin plate. Plate rest of media on 3 large (150x15mm) LB ampicillin plates per electroporation reaction. Grow at 37°C overnight.
7. Estimate the total number of colonies on all plates based on the 1:10,000 dilution plate. A rule of thumb is that you should target 100 times as many colonies as unique sequences that are in the oligo pool. If yield is satisfactory, recover colonies by adding 5 mL of LB to plates, and gently scrape surface of plates until colonies are no longer adhered to plate. Tilt the plate and recover the media with a pipette from the lower end of the plate. Repeat for each plate. Pellet cells by centrifugation at  ≥ 3,400 x g for 10 minutes and proceed with plasmid DNA purification using ZymoPURE™II Plasmid Maxiprep Kit. If yield is unsatisfactory, increase scale of previous steps and ensure DNA is free of contamination. This plasmid pool is referred to as “1°”.
8. PCR amplify intron 2 and exon 3 from modified pcAT7-Glo1 using Exon_3_MfeI_FWD and Exon_3_XbaI_REV for 35 cycles with 62°C annealing temperature. Purify PCR product with column-based purification. Digest purified PCR product with MfeI-HF and XbaI at 37°C for one hour and heat inactivate at 65°C for 20 minutes. Electrophorese DNA on 0.5% agarose gel and gel purify digested DNA fragment (~900 bp). If purified DNA has high contamination (260/230 < 2.0) after gel purification, use column purification or ethanol precipitation to further remove contaminants.
9. Digest at least 3 μg of 1° plasmid pool (product of step 7) with MfeI and SpeI at 37°C for one hour and heat inactivate at 65°C for 20 minutes. Add 1 μL CIP / μg of DNA and incubate at 37°C for 30 minutes. Electrophorese DNA on 0.5% agarose gel and gel purify digested DNA fragment (~5.8 kb). If purified DNA has high contamination (260/230 < 2.0) after gel purification, use column purification or ethanol precipitation to further remove contaminants.
10. Setup ligation in 1:3 vector to insert ratio with 500 ng of digested vector (product of step 9) and insert (product of step 8) in 85 μL total volume and incubate at 65°C for 5 minutes and put reaction on ice. Add 10 μL of 10x T4 ligation buffer and 5 μL of T4 ligase. Incubate ligation overnight at 16°C. Heat inactivate at 65°C for 20 minutes. Purify ligation using column purification and elute in 10 μL.
11. Transform and recover final plasmid DNA by repeating steps 6 and 7 with plasmid from step 10 to. This plasmid pool is referred to as “2°”.

## Plasmid Sequencing

1. Dilute a small volume of plasmid from plasmid assembly steps 11 and 7 to a final concentration of 10 ng/μL. Amplify diluted 1° plasmid with 1o_R1_FWD and Plasmid_R1_REV and amplify diluted 2° with 2o_R1_FWD and Plasmid_R1_REV for 13 cycles with annealing temperatures of 65°C and 63°C respectively.
2. Electrophorese products on 2% agarose gel and gel purify amplicons (~300 bp for 1° and ~360 bp for 2°).
3. PCR amplify R1 amplicons with R2_i5_FWD and R2_i7_REV and 10 ng template (product of step 2) for 10 cycles with annealing temperature of 72°C.
4. Electrophorese products on 2% agarose gel and gel purify amplicons (~370 bp for 1° and ~430 bp for 2°).

## RNA Sequencing

1. After harvesting RNA (as described in methods section), synthesize cDNA using up to 5 μg RNA per sample using Superscript III and VS_RT_UMI_primer using manufacturer’s protocol for a gene-specific primer.
2. PCR amplify cDNA using VS_FWD (one of VS_FWD1, VS_FWD2, or VS_FWD3) and R2_i7_REV for 10 cycles with an annealing temperature of 68°C and 5 cycles with an annealing temperature of 72°C. Note that the VS_FWD primers bind to the same sequence, but are offset by 0-2 nucleotides to add more diversity to the first few cycles of sequencing, which will result in higher sequencing quality.
3. Electrophorese products on 2% agarose gel and gel purify amplicons (between 300 and 500 bp)
4. Ensure no undesired small amplicons are present in samples using an Agilent TapeStation or equivalent equipment. If so, repeat gel extraction, or size select with AMPure beads. 
5. Combine samples from step 4 and plasmid sequencing samples at desired ratios for sequencing and prepare for sequencing on an Illumina platform with read 1 length of 75 nt and read 2 length of 225 nt. This configuration ensures that read 1 identifies the first junction, and read 2 identifies the second junction, the barcode, and the UMI.

## Sequences
| Primer | Sequence | Notes |  
|---|---|---|
| Oligo_FWD|CTGACTCTCTCTGCCTC|
| Oligo_REV|TAAACGGGCCCTCTAGA|
| Exon_3_MfeI_FWD|GTGTGGAAGTCTCAGGATCG|
| Exon_3_XbaI_REV|AACGGGCCCTCTAGAGC|
| 1o_R1_FWD|acactctttccctacacgacgctcttccgatctCCACTGACTCTCTCTGCCTC|Lowercase bind R2 primers
| Plasmid_R1_REV|gtgactggagttcagacgtgtgctcttccgatctAGCGGGTTTAAACGGGCCCT|Lowercase bind R2 primers
| 2o_R1_FWD|acactctttccctacacgacgctcttccgatctAGCAGCTACAATCCAGCTACCA|Lowercase bind R2 primers
| R2_i5_FWD|aatgatacggcgaccaccgagatctacac**_index_**acactctttccctacacgacgctcttccgatct|Replace **_index_** with 8 nt i5 sequencing index
| R2_i7_REV|caagcagaagacggcatacgagat**_index_**gtgactggagttcagacgtgtgctcttccgatct|Replace **_index_** with 8 nt i7 sequencing index
| VS_RT_UMI_primer|GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNNNNNGCTGATCAGCGGGTTTAAACG|
| VS_FWD1|AATGATACGGCGACCACCGAGATCTACACCGACTTTGACACTCTTTCCCTACACGACGCTCTTCCGATCTGGCAAGGTGAACGTGGATGAAG|Interchangeable with VS_FWD2 and VS_FWD3; use different ones between samples to increase sequencing diversity
| VS_FWD2|AATGATACGGCGACCACCGAGATCTACACAGGCTGTCACACTCTTTCCCTACACGACGCTCTTCCGATCTNGGCAAGGTGAACGTGGATGAAG|See VS_FWD1 note
| VS_FWD3|AATGATACGGCGACCACCGAGATCTACACATGCCGAGACACTCTTTCCCTACACGACGCTCTTCCGATCTNNGGCAAGGTGAACGTGGATGAAG|See VS_FWD1 note

