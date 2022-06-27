# SCIPhIN

## Run SCIPhIN

SCIPhIN expects the sequencing information to be passed in form of the well known mpileup format (http://www.htslib.org/doc/samtools.html). In order to generate such a file you need to align your fastq files to a reference and post process the result (e.g., following the instuctions here: https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery). 

In order to see all available options type

`sciphin -h`

Executing

`sciphin -o result --in cellNames.txt --seed 42 --im example.mpileup`

will run SCIPhIN using the cell names provided in *cellNames.txt* (same order as in the mpileup file). Note that *cellNames.txt* is a tab delimited file with the cell name in the first column and a cell type identifier in the second column. The cell type can be either *CT* (tumor cell), *CN* (control normal cell), or *BN* (control bulk normal). Note that SCIPhIN assumes a pileup against a reference and ignores positions with 'N' as reference.

## Examples

### Standard

The following command will run SCIPhIN with default settings and tell it to learn the zygousity (--lz), allow the loss of an allel (--ll) and allow parallel mutations (--lp). In addition, here we only run the MCMC for 10000 iterations. As input we need the input_spec.txt, telling sciphin what the input samples are, and the mpileup file of all samples.

`sciphin --lz 1 --ll 1 --lp 1 -l 10000 -o sciphin --in input_spec.txt --im example.mpileup`

After SCIPhIN finished there will be a sciphi.vcf with the mutation calls in VCF format. In addition, there will be a file called sciphin.probs, listing the probabilities for the different mutations type for each cell and mutation. Further a sciphin.gv file contains the most likely tree in dot format. With `dot -Tpdf sciphin.gc > sciphin.pdf` you can easily generate a PDF with the most likely tree representation.

### Hill-Climbing

Due to the MCMC tree estimation SCIPhIN may take a long time. We therefore offer the possibility to use a hill-climbing apprach rather than the MCMC scheme. In order to atcivate the hill-climbing mote and retrieve a point estimate you can use the option --mlm:


`sciphin --lz 1 --ll 1 --lp 1 --mlm -l 10000 -o sciphin --in input_spec.txt --im example.mpileup`

### Learning CHI

When infering the probabilities of loosing a mutation or observing the presence of a parallel muatation, noise is a major problem. In fact, random noise would not cause prolems, but in many cases the noise is not totally random, but instead, depends on the sequence content, for example. Therefore, one can specify a penalty for the loss of an allele (--llp) and the occurence of a parallel mutations (--lpp).


