# mt_contam_domestic_green
Script to estimate mitochondrial contamination in NGS data from ancient domestics following the Green 2008 approach (https://doi.org/10.1016/j.cell.2008.06.021)

This method first identifies sites in a dataset of diverse mitogenomes that are nearly fixed. Consequently, these sites, if different in a subject, can be used to estimate contamination from other samples. Non-consensus reads at such sites are divided by the total number of reads at the position to obtain a point estimate of mitochondrial contamination. We restrict this analysis to sites with at least 10x coverage. Furthermore, transition sites with a C or G in the consensus mitogenome are excluded to avoid over-estimation due to post-mortem damage. Standard errors are estimated assuming a binomial distribution around the point estimate.

Note that the mitochondrial contamination can be very different from nuclear contamination (e.g. due to different numbers of mitochondrial between cell tyes) and that NUMTs can be problematic due to spurious mapping to the mitogenome.

## Preparation of MT dataset

*We outline how we used the scripts here but you can obviously collect your own dataset of mitogenomes.*

We first downloaded the sequences used for Dometree (https://doi.org/10.1111/1755-0998.12386) from DRYAD: https://datadryad.org/stash/dataset/doi:10.5061/dryad.cc5kn and unziped the file.

The file contains folders for cattle, dogs, goats, horses, pigs, sheep, yaks and chickens. To collect sequences for our species of interest, we went into the corresponding folder and extracted the second sequence from each pairwise alignment:

``for file in $(ls *.fasta); do tail -n 2 $file >> all.fasta; done``

We then remove gaps from the file:

``sed -i "s/-//g" all.fasta``

Before performing a multiple sequence alignment it is important to add a reference sequence, the MT sequence that has been used for mapping short reads, to the FASTA file. We did this manually since the Dometree database had used a different reference genome.

Next, one needs to perform a multiple sequence alignment with the aligner of choice, e.g. MAFFT, Muscle or Clustal. Produce the output as a FASTA file.

**Make sure that the resulting alignment has the reference sequence at the first position, the following script assumes that. Otherwise the identification of potentially private sites will be messed up. Some aligners change the order of sequences compared to the input file.**

## Preparation of site list

Now, we need to calculate the minor allele frequency (MAF) at each position relative to the reference sequence. We use the Python script `green_mt2tab.py` for this step. It requires to specify the input file (the multiple sequence alignment FASTA) with `-i` one can also provide custom MAF file names with `-o`.

This step needs to be conducted only once per MSA of mitogenomes, the MAT/sites file can then be re-used to estimate contamination in multiple samples.

## Estimating conamination per sample

Next, we can estimate contamination for our sample. We need a BAM file, the reference genome FASTA and `samtools`:

``samtools mpileup -B -q 30 -Q 30 -r MT -f $ref_seq $bamname | ./mt_contamination.py -s $MAF_file -f 0.05``

`-s` specifies the name of the sites/MAF file. With `-f` the MAF cutoff can be changed, the default is 0.05. The following results are printed to `stdout`: point estimate, informative sites, consensus alleles, total alleles, lower end CI, upper end CI
