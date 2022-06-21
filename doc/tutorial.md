# *Candida albicans* SDG assembly 
This is a short walkthrough to illustrate how SDG can be used to generate haplotype-specific chromosome-scale assemblies from Illumina short reads and Nanopore long reads. The example dataset is *Candida albicans* from the following publication;

High-Quality Genome Reconstruction of Candida albicans CHN1 Using Nanopore and Illumina Sequencing and Hybrid Assembly.
Authors: Shipra Garg, Piyush Ranjan, John R. Erb-Downward and Gary B. Huffnagle.
DOI: [https://doi.org/10.1128/MRA.00299-21](https://doi.org/10.1128/MRA.00299-21) 

## 1. Download SDG
Download and compile the latest version of SDG from [https://github.com/bioinfologics/sdg](https://github.com/bioinfologics/sdg)

## 2. Download read datasets and check kmer histogram
The raw reads are available under accession numbers [SRX9854709](https://www.ncbi.nlm.nih.gov/sra/SRX9854709) (Illumina) and [SRX9854710](https://www.ncbi.nlm.nih.gov/sra/SRX9854710) (Nanopore) within BioProject [PRJNA692229](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA692229).

Once you have the FASTQ files, use the [K-mer Analysis Toolkit (KAT)](https://github.com/TGAC/KAT) to generate a kmer histogram from the paired-end reads. This will provide us with the unique coverage frequency estimate we need.
```
kat hist -o chn1_pe -t 8 -m 27 -H 100000000 SRR13441295.1_?.fastq
```

![](https://i.imgur.com/4wcK60X.png)

The histogram shows a fundamental peak at a frequency of 67 which contains unique, heterozygous content from the reads. This is followed by a further peak at a harmonic frequency of 134, containing homozygous content. This is a typical diploid spectra signature with high homozygosity. 

## 3. Create SDG datastores from the reads
You need to convert the FASTQ read files to datastores, one for the Illumina reads and one for the Nanopore reads. This is the format that SDG requires.
```
sdg-datastore make -1 SRR13441295.1_1.fastq -2 SRR13441295.1_1.fastq -t paired -o cachn1_pe
sdg-datastore make -L SRR13441294.1.fastq -t long -o cachn1_nano
```

## 4. Build contig graph
### The importance of the K-mer Compression Index (KCI) in graph construction
For a given node in a graph, the KCI indicates how many times that node should appear in the final assembly. It is calculated by taking unique kmers from the node and identifying how many times these kmers appear in the reads. In a diploid, the coverage of homozygous components will be twice that of the heterozygous components. A KCI close to 1 indicates a heterozygous node and a homozygous node will have a KCI close to 2. At higher ploidy levels, the KCI will generally indicate the number of haplotypes the node represents (eg. homozygous across 3, 4, 5, etc. haplotypes) provided the heterozygous coverage can be ascertained. In a triploid, a homozygous locus will have a KCI close to 3, but heterozygous regions will appear as a mix of 2:1 paired homologs, or 1:1:1 homologs. As the KCI is calculated genome-wide, it can indicate ploidy levels across the genome and be employed to check homozygosity, heterozygosity, ploidy, and erroneous assembly at a local scale. This makes manual analysis of the structure of the assembly graph possible. Haplotypes can be traced through the graph and ambiguous links between nodes manually checked for KCI support. Misassemblies caused by the crossing over of linked nodes from one haplotype to another can be identified and further corroborated by long-read support.

### Step 1: Create DBG from short reads
The SDG pipeline uses the KCI to create a de-bruin graph suitable for haplotype reconstruction in five major steps. In the first step, the initial graph is constructed from the reads. The graph topology is categorised into three groups. Tips are path dead-ends, which contain only an input connection from a previous node. Canonical repeats are repetitive nodes with a minimum of 2 input and output connections. These may be higher coverage areas that are found on more than one of the haplotypes or repetitive content. Finally, bubble sides are parallel nodes connecting to the same input and output nodes with no other connections, forming distinctive bubbles in the graph. Depending on your dataset and available resources this step can take quite a while. The test dataset should run in around 4 hours.
```
01-contigger/01-dbg.py -o cachn1 -c 3 --kci_k 51 -p cachn1/cachn1_pe.prseq
```
Note: This step takes the longest to compute, up to several hours.

### Step 2: Remove low coverage nodes
In the second step, unique anchors are identified as nodes which occur uniquely in the genome and are represented by a KCI of 1. Low KCI nodes and tips are removed if the nodes they are connected to have an alternative route through the graph. Repeats are retained for phasing, but some erroneous nodes will remain.
```
01-contigger/02-clean.py -o cachn1 -u 67 --low_coverage 40
```
### Step 3: Expand short repeats
The third step expands the canonical repeats by choosing paths with the highest support. This connects bubble sides into parallel phased nodes through the canonical repeat. The paths with low support from a user defined noise threshold can be disregarded.
```
01-contigger/03-short_repeats.py -o cachn1 -u 67 --pe_rep_max_noise 5
```
### Step 4: Run Strider for further repeat resolution
The fourth step introduces STRIDER with long-read support to resolve longer contigs. A path is mapped between the unique bubble side nodes and the repeats, expanding out these regions into haplotype specific nodes. As STRIDER is run consecutively, these regions become further expanded as more of the graph is solved. 
```
01-contigger/04-strider.py -o cachn1 -u 67
```
### Step 5: Canonical repeat resolution using long reads
The final step uses the long-reads to give a final level of haplotype expansion. This traces more paths through the canonical repeats and gives support to paths between canonical repeats and bubble side nodes. 
```
01-contigger/05-long_repeats.py -o cachn1 -u 67 -l cachn1/cachn1_nano.loseq --lr_min_support 5 --lr_max_noise 5 --lr_snr 5
```
## 5. ReadThreadGraph
SDG's haplotype-specific assemblies are constructed from a ReadThreadsGraph. This graph contains single-copy anchor sequences, extracted from an assembly graph; and threads representing long read alignments to these anchor sequences. SDG then produces a haplotype-specific partition of the graph, where each disjoint class contains a subset of its anchors and threads representing a single haplotype in a unique region. The anchors of each class can be then ordered, producing a backbone that is used to assemble the haplotype-specific consensus.

```
02-scaffolder/06-split_and_map.py -i cachn1_05_long_repeats.sdgws -o cachn1 -u 67 -s 1500
```
## 6. KaryoThreader
Using long read mappings to the anchor graph a ReadThreadsGraph is constructed, with threads representing each long read as a linkage of anchors. Clusters of anchors and threads are created, representing haplotype-specific regions of the genome as anchors and the threads that join them. The TotalSorter algorithm subtracts each stable partial solution from the problem, aiding in progressively cleaning up the remainder of the problem until the whole set of anchors and threads is divided into partial solutions.

Properties of the ReadThreadsGraph can be used to cleanup and define proximity between pairs of nodes, pairs of threads, etc.

Two fundamental functions are the Thread neighborhood (i.e. the group of trheads that share nodes with a thread) and Node neighborhoouds (i.e. the group of nodes that share nodes with a node).

Thread neighborhood can be computed by simply going thorugh the nodes and incrementing a counter for each pair of threads seen concurrently.

For node neighbourhood, as to simplify computation, the neighborhood is restricted to nodes seen around a set number of nodes around in a thread. I.e. for the node in position x in a thread, all nodes in positions x-radius and x+radius have their neighborhood values incremented.

Long read mapping was used to scaffold haplotypes. The SDG HappySorter uses long-reads to bridge unlinked nodes, mapping only to the unique anchors of the long contigs. This final level of assembly can provide support for complex regions that are otherwise difficult to resolve, such as highly repetitive regions and regions with low coverage. It uses the unique nodes from STRIDER to map against the long reads, supporting further haplotype resolution.

```
02-scaffolder/07-thread_and_scaff.py -o cachn1 -u 67 -s 275 --min_hits 3 --min_kci .25 --max_kci 2.75 --min_links 10 --max_thread_count 500 
```
## 7. Manual analysis
From this point, the graph can be accessed for manual analysis. This step resolves the final complexities that escaped the parameterisations for reasons such as low/high coverage, repetitive, misjoins, etc. The contig graph can be used to access all nodes and linkages, while the reduced graph shows us the longer-scale structure. 


[Bandage](https://rrwick.github.io/Bandage/) can be used to explore the graph topology.

Nodes and links can be included or excluded by creating a file of inclusions/exclusions and using the flag --rrtg_edits. For example:

```
exclude node 2253693
include node 31573
include link 2165929 2148120
exclude link 2165456 2223485
```
```
02-scaffolder/07-thread_and_scaff.py -o cachn1 -u 67 -s 275 --min_hits 3 --min_kci .25 --max_kci 2.75 --min_links 10 --max_thread_count 500 --rrtg_edits cachn1.rrtg_edits
```

![](https://i.imgur.com/pCSmfKY.png)
A canonical repeat (red) causing a misjoin between two chromosomes. This node can be manually removed, which will uncouple the chromosomes. 

## 8. Join and fill
The join and fill function fills in gaps in the scaffold graph with overlapping nodes, removing transitive links and creating contiguous sequences. This is used to back fill as much content as possible, leaving out of the assembly only repeat regions that were not resolved.

![](https://i.imgur.com/p3FN0Qh.png)
All 8 chromosomes fully resolved and mapped back to the published assembly from the example dataset, so that each chromosome has an indivdual colour. Chromosome R still has a missing linkage to a small cluster of nodes that could not be reconnected under the parameters used. 
