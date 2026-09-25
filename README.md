# DEVINATE

## Long-reads TEs insertions finder 

### Installation<a name="Installation"></a>

#### Requirements<a name="Requirements"></a>

For this very first version of DEVINATE, we strongly recommend using the provided Conda environment to install all required dependencies :

```bash
conda env create -f devinate.yml
```

Then activate it:

```bash
conda activate devinate
```

| Software   | Version |
|------------|---------|
| Python     | 3.12    |
| Biopython  | 1.87    |
| pysam      | 0.23.3  |
| SeqKit     | 2.9.0   |
| minimap2   | 2.27    |
| samtools   | 1.21    |
| BamDash    | 0.5.0   |

DEVINATE was designed on these tools versions. 

#### Using Git<a name="git"></a>

Once the requirements fullfilled, just *git* clone

```bash
git clone https://github.com/paillerv/DEVINATE.git
```

### How to use this script?

The main bash script DEVINATE_main.sh will launch sequentially the python sub-scripts.

```bash
bash DEVINATE_main.sh \
   ./reads.fasta.gz \
   TE-consensus.fasta \
   reference-genome.fasta \
   OUTPUT_name
```
