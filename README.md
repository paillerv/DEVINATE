# DEVINATE

## Long-reads TEs insertions finder 

- [Requirements](#requirements)
- [Installation](#installation)
    - [Using Git](#git)
    - [Using Singularity](#singularity)
- [Configuration](#configuration)
- [Usage](#usage)

#### Requirements<a name="requirements"></a>

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

### Installation<a name="installation"></a>

#### Using Git<a name="git"></a>

Once the requirements fullfilled, just *git* clone

```bash
git clone https://github.com/paillerv/DEVINATE.git
```

DEVINATE was developed using these specific tools and versions. This may change in the near future, with a Singularity option. 

#### Using Singularity<a name="singularity"></a> 

TO DO. 

### Configuration<a name="configuration"></a> 




### Usage<a name="usage"></a> 

The main bash script DEVINATE_main.sh will launch sequentially the python sub-scripts.

```bash
bash DEVINATE_main.sh \
   ./reads.fasta.gz \
   TE-consensus.fasta \
   reference-genome.fasta \
   OUTPUT_name
```
