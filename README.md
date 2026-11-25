# WES and molecular dynamics analysis pipeline for the identification of the novel missense variant in the USP26 gene in asthenoteratozoospermia

Code accompanying:

> Mousavi SZ, Askari M, McElreavey K, Bashamboo A, Salehi N, Rastari M, Rokhsattalab Z, Mohammad Soltani B, Totonchi M.  
> A novel candidate missense variant in the catalytic domain of USP26 associated with asthenoteratozoospermia.  
> *Gene*. 2025.  
> https://doi.org/10.1016/j.gene.2025.149945

## What this repository contains

Computational analysis related to a missense variant in the catalytic domain of USP26. The repository has two parts:

| File | Content |
|---|---|
| [`wes_pipeline.sh`](wes_pipeline.sh) | Whole-exome sequencing (WES) variant-calling and annotation pipeline (FastQC, BWA-MEM, samtools, GATK 4.1.9.0, ANNOVAR; reference build hg19) |
| [`gromacs-minimization.md`](gromacs-minimization.md) | GROMACS commands for system setup, energy minimization and NVT/NPT equilibration of the USP26 protein |

## 1. WES pipeline (`wes_pipeline.sh`)

Steps in the script, in order:

1. Read quality control (FastQC)
2. Alignment to hg19 (BWA-MEM), SAM to BAM conversion and filtering of unmapped reads (samtools), coordinate sorting (GATK SortSam)
3. Read groups and duplicate marking (GATK AddOrReplaceReadGroups, MarkDuplicates)
4. Base quality score recalibration (GATK BaseRecalibrator, ApplyBQSR, AnalyzeCovariates)
5. Per-sample GVCF calling (HaplotypeCaller), then family-level joint calling (CombineGVCFs, GenotypeGVCFs)
6. Variant quality score recalibration for SNPs and indels (VariantRecalibrator, ApplyVQSR)
7. Annotation (ANNOVAR, hg19)

**Requirements.** bwa (BWA-MEM), samtools, FastQC, GATK 4.1.9.0, ANNOVAR, and an hg19 reference FASTA (index once with `bwa index`).

**Usage.** The script is a template. Edit the variables at the top (`SAMPLE`, `FAMILY`, `THREADS`, and all `.../` paths to the reference genome, GATK resource VCFs and ANNOVAR databases), set `FAMILY_GVCFS` to the per-sample GVCFs of the family, and run:

```bash
bash wes_pipeline.sh
```

Reference genome, GATK resource files and ANNOVAR databases are not included and must be downloaded separately.

The alignment step follows the original four-FASTQ layout: `ForwardPairs.fastq` / `ReversePairs.fastq` are aligned as a paired-end set, and `ForwardUnpairs.fastq` / `ReverseUnpairs.fastq` (leftover unpaired reads after trimming) are each aligned as single-end. `FAMILY_GVCFS` in the script contains repeated placeholder entries; replace them with one distinct GVCF per family member.

## 2. Molecular dynamics (`gromacs-minimization.md`)

The commands for the simple protein simulation (box setup, solvation, ion addition, energy minimization, NVT and NPT equilibration) are in [`gromacs-minimization.md`](gromacs-minimization.md).

- This is the **standard GROMACS workflow** as described in the public GROMACS lysozyme tutorial by J. A. Lemkul (<http://www.mdtutorials.com/gmx/lysozyme/index.html>). It is not a novel method.
- **Parameter files are not included.** The `.mdp` files used by the commands (`ions.mdp`, `minim.mdp`, `nvt.mdp`, `npt.mdp`) are not in this repository. Standard `.mdp` files can be obtained from the tutorial above; they must be adapted to the force field and water model you select in `pdb2gmx`.
- Input structure: the commands use `USP26_model.pdb` as a placeholder name for the USP26 structure file.
- The documented protocol ends after NPT equilibration; a production run is not documented here.

## Data policy

No identifiable patient data. Variant-level and public annotations only. Raw sequencing files, BAM/VCF files and other intermediates are excluded via `.gitignore`.

## License

MIT License, see [`LICENSE`](LICENSE).

## Citation

Please cite the paper above.
