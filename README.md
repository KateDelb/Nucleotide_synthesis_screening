## Identifying function-relevant signatures in protein models for biosecurity screening

I propose a (very) simple screening pipeline for a nucleic acid sequencing screening dealing with an unknown sequence:
"ATGAATCAACATAATACCCAAATCAATAAATTTATCTTTTTTAGTAGCTTCAGAATCATCAGTACCACCCCAATTCAATTTATCAACAATAATAGGACCAGCAACAGCAATAGCAACATATTCATCATTAGAAGCAGGAGAATTCAAAACAATACCAACAGCTTTTTTATTATCAGCAGTATTCCAAGTTTTACAAACAGCACCATTAGAATTACCTTTTTCACAATCACCAGCATGCAAATTATCAATACCAGAATTTTTATGAGATCTAATATGAACCAAACCCAT"

main.py describes the following steps:

-> Identify if this is eukaryote/viral/bacterial. If this is eukaryotic it is less likely to be highly dangerous.
NOTE: The classify_organism function is a bare minimum check of GC and CpG content and can not tell Eukaryote from Viral for the moment. 

This seems to be an ORF (starts with ATG), would be interesting to 
-> Check if there is a STOP codon somewhere. If there isn't, it *could* plausibly imply that other parts of the ORF are ordered at another company to avoid being flagged.
NOTE: check_orf function is simplified: assumes 1 start codon preceding 1 standard end codon (TAG, TAA or TGA)

-> Run BLASTN (I assume this will be negative if the initial internal screening didn't bring up anything)

-> Run BlASTx: See if we can identify an known homologous protein sequences, if this can be linked to any known/dangerous sequences

-> Predict protein structure to see if any functional class can be determined using HMMer (polymerase domain architecture, toxin-like folds?)


Obviously also: 
- Check who is ordering, anonymous or non-researcher, no clear credentials should immediately be flagged 
- what quantity is this being ordered?
- What is the justification of the research to be done

Things that could additionally be checked:
- update classify_organism function (low prio, not very informative)
- Check against curated threat databases (like GISAID for priority pathogens), check for structural viral proteins (fairly low risk in se but could add clues if already suspicious)


Disclaimer:
Claude Haiku was used for the code of classify_organism, blastn, blastx and hmmer with minor corrections as well as for the suggestion of the HMMer tool to identify known domains.