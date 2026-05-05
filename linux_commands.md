## linux commands - reproducible workflow 

we start our analysis with the aim of building a profile Hidden Markov Model for the Kunitz-type protease inhibitor domain by filtering the search on PDB thereby recruiting monodomains of Kunitz type according to the given constraints:

- sequence length above 40 aminoacids
- sequence length below 80 aminoacids
- PFAM identifier being PF00014
- Data collection resolution equal or below 3 Ångström

From the filtered search on the advanced query builder in PDB we obtain 135 entries; before proceeding with the download of the research's output we create a costum Tabular Report (available in this repository) highlighing the sequence identifier(entry ID), the sequence and the chain.

Although the advanced search targeted the PF00014 identifier, structural data often contains multiple chains per entry. For homodimeric structures, where Chain A and B are sequence-identical, only a single representative chain was retained to avoid sequence bias. For heterodimeric complexes involving auxiliary fragments (e.g., Chains Q or J), a length-dependent filter was applied to select only the primary Kunitz monodomain. This curation resulted in a refined dataset of unique, high-resolution sequences optimized for HMM construction.

---
### data cleaning step with bash commands

By downloading the Tabular report (provided) in the CSV format,  the following commands can be reproduced in order to obtain such clean search:

cleaning the CSV by removing lines starting with commas, stripping quotes, and converting the format into a space-separated list

```
grep -v "^," custom_tabular_report.csv | tr -d '"' | tr "," " "
```

To ensure the Profile HMM was built on high-quality, non-redundant monodomains, a manual curation step was performed. the output was redirected to the pdb_seqs.txt
```
awk -F "," '{if ($1!="" && length($2)<80 && length($2)>40 ) print $1,$3,$2}' rcsb_pdb_custom_report_20260415204951.csv | tr -d '"' | awk '{print ">"$1$2;print$3}' > pdb_seqs.txt
```

Then the output was checked to verify the amount of lines obtained 
```
grep ">" pdb_seqs.txt | wc
```
---

### removing the redundancy and clustering with MMseqs2

MMseqs2 toolkit removes the redundancy by using a local alignment with the parameters being for our analysis:

- minimum sequence identity 0.95
- minimum alignment coverage 0.90

The output is a cluster of sequences which for our input given(pdb_seqs.txt) retrieved 19 Clusters 
We proceed by creating a file.clust in which we store the output of MMseqs2

```
vi pdb_seqs.clust
```
and we copy and paste inside it our values; then we check for the presence of each cluster

```
grep -A 1 "^Cluster#" pdb_seqs.clust | grep ">" | wc
```

we keep only a representative sequence for each cluster (identified as the first per cluster) redirecting our output to representative ids togheter with the chain

```
grep -A 1 "^Cluster#" pdb_seqs.clust | grep ">" | tr -d ">" > pdb_id_chain.rep
```

since for our structural alignment we will not need the chains we clean our pdb_id_chain.rep

```
cut -c 1-4 pdb_id_chain.rep > pdb_id.rep
```
---

### generating the structural alignment 

There are different methods to perform the structural alignment such as using *PDBeFold* on your machine or run it web. The choice made was tu use the *mTM-align* Multiple Protein alignment tool by installing it on the local machine. All of the instructions needed to correctly install the mTM-align can be found inside the 
```
cd mTM-align  / ls / cat README.txt 
```
execute and install

In order to perform the structural alignment starting from our representative PDB ids we have to seek the structures in PDB. This step can be perfomed using the following for loop ; the structures of interest can be either store in our local directory or create a dedicated directory with the *mkdir* command for better organization 
```
mkdir pdb_chains
```

```
for i in $(cat pdb_id.rep); do wget "http://files.rcsb.org/view/${i}.pdb"; done
```

Once the PDB files are downloaded, we must isolate the specific structural coordinates for the Kunitz domain. In a PDB file, the relevant structural data is contained in lines starting with the ATOM record. Within these lines, position 22 (column-based) identifies the Chain ID.

To automate the extraction of specific chains from each downloaded PDB file, we follow this workflow starting from our original file * pdb_id_chain.rep* containig our chain of interest:

```
sed 's/\(.\)$/\t\1/' pdb_id_chain.rep > pdb_ids_tab_separated.rep
```

We used the Vim editor to develop a Python utility designed to parse the raw PDB files and generate new, filtered files containing only the relevant structural data.

```
vi get_chain.py
```

The file get_chain.py can be found in this repository 

```
while IFS=$'\t' read -r pdb chain;do python get_chain.py $pdb".pdb" $chain > ""pdb_chains/$pdb$chain".pdb"; done < pdb_ids_tab_separated.rep
```

in order to perform the structural alignment mTM-align wants the complete path inside one unique file so we perform

```
ls pdb_chains/* > my_input_mTM.txt
```

we run the alignment using mTMalign

```
mTM-align/src/mTM-align -i my_input_mTM.txt > my_output_mTM.txt
```

we check visually that our output has aligned all of our structures (19 in my case) and we store the output in a new file *kunitz.ali* adding the header symbol of the  owing to FASTA files formatting 

```
tail -n +26 my_output_mTM.txt | head -n 19 | awk '{print ">"$1;print $2}' > kunitz.ali
```
inspecting the output 
```
less kunitz.ali
```

After an initial multiple sequence alignment (MSA) and visual inspection of the output, a final manual curation step was performed to ensure the highest statistical confidence for the HMM profile.

Two specific entries were identified as outliers:
- 1D0D (Chain A)
- 5YW1 (Chain A)

Despite meeting the initial length (40–80 AA) and PFAM (PF00014) criteria, these sequences showed significant divergence. Most importantly, they did not respect the Kunitz domain constraint of the six conserved cysteines, which are responsible for the three essential disulfide bridges stabilizing the fold. To obtain the cleanest possible signal and a profile with higher biological relevance, these two sequences were removed from the input dataset.
```
cp my_input_mTM.txt my_clean_input_mTM.txt

rm 1D0DA.pdb 5YW1A.pdb my_clean_input_mTM.txt

mTM-align/src/mTM-align -i my_clean_input_mTM.txt > my_clean_output_mTM.txt
        
tail -n +26 my_clean_output_mTM.txt | head -n 19 | awk '{print ">"$1;print $2}' > kunitz_clean.ali
```

---

### generating the hidden markov model 

After having performed the structural alignment we have to transform this in a sequence alignment in order to understand possible issues that are present in the structural alignment.

```
hmmbuild kunitz_clean.hmm kunitz_clean.ali
```

by inspecting the output of the alignment we see that matches start at  position 10 and end at position  69, indicating an average length of 57 aminoacids. This aligns perfectly with the biological reality of the PF00014 family, which typically spans between 50 and 60 amino acids. This means that the alignment is now very "tight." Only 4 positions are being treated as insertions/deletions, which indicates a very high-quality Multiple Sequence Alignment (MSA).

```
less kunitz_clean.hmm
```

---

### Validation of the model against uniprot

in Uniprot we perform an advanced search to create a reference dataset that we can use to benchmarking our method; we divided therefore our search between:
- reviewed true and pfam PF00014
- reviewd true not pfam PF00014

download the output in .FASTA format 
*i check they have the same entries they have in Uniprot*
```
grep ">" uniprotkb_reviewed_true_NOT_pfam_PF00014.fasta|wc
grep ">" uniprotkb_reviewed_true_PF00014.fasta|wc
```
We run the search on the two downloaded outputs 

#### checking for the True positive (TP)
```
hmmsearch --max --noali --tblout positive_kunitz.search -Z 1000 kunitz.hmm uniprotkb_reviewed_true_PF00014.fasta
```

#### checking for the False positives (FP) - meaning thought to be non-kunitz but are actual kunitz
on the clean search 
- there are two e-values (1. full sequence/ 2.best domain)
- there is also the option in hmmsearch of --max to remove the heuristic and get an optimal alignment
```
hmmsearch --max --noali --tblout negative_kunitz.search -Z 1000 kunitz_clean.hmm uniprotkb_reviewed_true_NOT_pfam_PF00014.fasta
```


To test the performance and make the optimization we run a cross validation procedure; the easiest case is to separate the dataset in two subsets
- you can do this initially and then run the calculation later 
- or you can split the dataset and run hmmsearch on the two datasets
- fix the option in hmmsearch that fixes the number of queries you are executing - hmmsearch { hmmsearch -h}
 -Z <x>        : set # of comparisons done, for E-value calculation

to be on the safe side you can tell hmmsearch to run the search on the two databases a fixed amount of times so that even if the two dataset are of different sizes it is still possible to run the same amount of tests

#### creating the list of matches 

*to get only the identifier and the e-value we want to parse the output of the HMMsearch*

To this respect here we add a fixed column in third position for binary classification which is a format compatible with Performance Evaluation Metrics.

- The "1" (Positive Label): Indicates the sequence belongs to the "True Kunitz" (reviewed PF00014) reference set.

- The "0" (Negative Label): Indicates the sequence belongs to the "Non-Kunitz" (reviewed not PF00014) control set.

**for the negative search**
we first redirect the identifier - evalue and binary laber to the match file
```
grep -v "^#" negative_kunitz.search | awk '{print $1"\t"$8"\t0"}' > negative_kunitz.match

```
we extract only the identifier for the matches 
```
awk '{print $1}' negative_kunitz.match | sort > negative_kunitz_match.ids
```
we extract all initial database identifiers to know the original population

```
grep ">" uniprotkb_reviewed_true_NOT_pfam_PF00014.fasta  | awk '{print $1}'| tr -d ">" | sort > negative_kunitz.ids
```

we use the comm command to find IDs present in the database so the initial negative_kunitz.ids but not in the search result 

```
comm -23 negative_kunitz.ids negative_kunitz_match.ids | awk '{print $1"\t100\t0"}' > negative_kunitz.nomatch
```

note: We assign an E-value of 100 to these sequences. In a HMMER context and particularly in our case after having investigated the maximum e-values, an E-value of 100 is statistically insignificant, effectively telling your downstream analysis scripts that the model "rejected" these sequences.

Finally we merge the hits and the non-hits so that the final file size matches your original UniProt search.

```
cat negative_kunitz.match negative_kunitz.nomatch | sort -R > total_negative_kun.matches
```

**for the positive search**

```
grep -v "^#" positive_kunitz.search | awk '{print $1"\t"$8"\t1"}' > positive_kunitz.match
```

To strictly evaluate the HMM's performance on unseen data, we must remove sequences used during the model construction (the training set) from the positive validation dataset.

We first extract the PDB accessions used in the initial structural alignment and convert them to a standardized uppercase format:

```
cut -d'/' -f2 my_clean_input_mTM.txt | cut -c1-4 | tr '[:lower:]' '[:upper:]' | sort -u > training_pdb_codes.txt
```

Since the training set consists of PDB structures and the validation set consists of UniProt entries, we perform a sequence-identity search to find exact matches.

* ID Mapping: PDB codes(training_pdb_codes.txt) are mapped to UniProt sequences (downloaded as id_mapping.fasta)

* Database Creation: A ```local BLAST database``` is generated from the reviewed UniProt Kunitz dataset.

* Sequence Alignment: We run ```blastp``` to identify sequences in our validation set that are identical to our training set.

We isolate only those hits with 100% identity and an E-value of 0.0, ensuring we only remove exact duplicates:
```
awk '$3 == 100.000 && $11 == 0.0 {print $2}' results_blast.txt | cut -d'|' -f2 | sort | uniq > ids_to_rm.txt
```

The final "positive" match file is cleaned by removing these identifiers, resulting in a purified test set (e.g., 394 sequences) ready for unbiased cross-validation:
```
grep -v -f ids_to_rm.txt positive_kunitz.match > positive_kunitz_clean.match
```
---

### k-fold cross validation

Dataset partitioning 

The purified positive set (n=394) and the reconstituted negative set (n=574,229) were split into two equal subsets (Set 1 and Set 2). To ensure a balanced evaluation, the negatives (including the imputed "non-hits" with arbitrary E-values of 100) were shuffled before splitting to prevent any order-based bias.

```
head -n 197 positive_kunitz_clean.match > kunitz_set_1.txt
tail -n 197 positive_kunitz_clean.match > kunitz_set_2.txt

head -n 287115 total_negative_kun.matches >> kunitz_set_1.txt
tail -n 287114 total_negative_kun.matches >> kunitz_set_2.txt
```

We iterated through a range of E-value thresholds using the *performance.py* script. The ```Matthews Correlation Coefficient (MCC)``` was used as the primary metric for optimization because it is the most robust measure for highly imbalanced datasets (where negatives far outnumber positives).

```
for i in `seq 1 15`; do python performance.py kunitz_set_1.txt 1e-$i; done > kunitz_set_1.performance

for i in `seq 1 15`; do python performance.py kunitz_set_2.txt 1e-$i; done > kunitz_set_2.performance
```
Ultimately we sort our values according to the Matthew Correlation coefficient to understand the best e-value threshold for our model 

```
sort -k6,6gr kunitz_set_1.performance
sort -k6,6gr kunitz_set_2.performance
```
