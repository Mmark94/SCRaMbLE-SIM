# SCRaMbLE-SIM: Genome SCRaMbLE evolution simulator
<!--
[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/wouter_decoster.svg?style=social&label=Follow%20%40Mm94Marco)](https://twitter.com/Mm94Marco)
-->

## Table of Contents
- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Applications](#applications)
- [Usage](#usage)
    - [Genome evolution through SCRaMbLE](#genome_evolution)
    - [Simulation of long sequencing reads](#reads_SIM)
    - [Genome composition and segment copy number analysis](#CN)
    - [SCRaMbLE population simulator](#population)
    - [Mortality rate analysis of SCRaMbLE in a populatio](#mortality)
    - [Chromosome size evolution analysis](#chr_size)
- [Getting help](#help)
- [Citation](#citation)
- [Limitations](#limitations)

## <a name="introduction"></a>Introduction
SCRaMbLE-SIM is a genome simulator software that simulates genome evolution through genomic rearrangements like deletions, inversions, duplications and translocations.
This project is inspired by the Sc2.0 technique SCRaMbLE which allows stochastic rearrangements between symmetrical loxP sites within the synthetic chromosomes.

## <a name="dependencies"></a>Dependencies
- Biopython

## <a name="applications"></a>Applications
SCRaMbLE-SIM can be applied to answer both practical and theoretical questions in the fields of Genome Evolution and Bioinformatics.
The areas where we think this software could be most useful are: modelling Genome Evolution and Genome Size Evolution, benchmarking software and algorithms for Genome Evolution, genetic distance and Phylogenesis, for benchmarking and validating software and algorithms for *de novo* assembly, mapping, and structural variant calling.

These are some applications that were tested so far in this project and their usage is described in the following section.
- Synthetic genome evolution through SCRaMbLE
- Genome evolution of any genome through genomic rearrangements
- Simulation of long sequencing reads
- Genome composition and segment copy number analysis
- SCRaMbLE population structure analysis
- Mortality rate analysis of SCRaMbLE in a population
- Chromosome size evolution analysis


## <a name="usage"></a>Usage
Following are some applications of SCRaMbLE-SIM with some examples of how to use them.

### <a name="genome_evolution"></a>Genome evolution through SCRaMbLE
All the essential functions to SCRaMbLE a chromosome are in the scripts SCRaMbLE_simulation_3.py or SCRaMbLE_simulation_3_circular.py. You can SCRaMbLE a chromosome within python by calling the function **force_SCRaMLE_lin_cir** or **SCRaMbLE_muliple_chrs** if you have multiple chromosomes.
This function takes as input the initial chromosome path sequence (the list of segments. For example, [1,2,3,...,44]), and the number of SCRaMbLE events to simulate. SCRaMbLE-SIM can simulate SCRaMbLE in one or many linear or circular chromosomes.

Moreover, each rearrangement type can have a different probability of happening. The default probability settings are respectively [NULL, DELs, INVs, DUPs]=[0, 2, 2, 1].

SCRaMbLE-SIM also allows the user to determine the length distribution of the SCRaMbLE events. Based on evidence from our previous works (Z. Luo et al., 2021; Y. Shen et al., 2016), the probability of a SCRaMbLE event decreases with the distance between loxPsym sites; therefore, smaller events are more likely than larger ones. To simulate this event distribution, we have chosen the event length randomly using a "Discretized Half-normal distribution". The default standard deviation 𝜎 is 10 LUs (10 LUs ≈ 29 Kb).

In addition, when multiple chromosomes are present, the probability of a SCRaMbLE event in one chromosome is proportional to the chromosome length. Moreover, in the case of multiple chromosomes, also translocations are possible and the probability of translocation can be set by the user and by default is 5 %.

SCRaMbLE-SIM allows the simulation of SCRaMbLE with essential LUs and the centromere LU (which must be present in any chromosome exactly once). All essential LUs must be present in the genome in at least one copy, if they are deleted, the program discards the event and generates a new one.

```py
force_SCRaMLE_lin_cir_events(syn_chr: list, Number_events: int, essential=[], circular=False, mu=0, sigma=7, CEN=[], force=True, probability=[0, 2, 2, 1], event_type=False)
```

I am also working on the script SCRaMbLE_DNA_chromosomes.py which will take as input a genome in fasta file and output a SCRaMbLEd genome in a fasta file.


### <a name="reads_SIM"></a>Simulation of long sequencing reads
SCRaMbLE-SIM can also simulate long sequencing reads such as those from Nanopore or PacBio by cutting input chromosomes. The main advantages of simulating reads are the complete control over the coverage (how many reads and LUs) and the read length distribution of the simulated reads.

To simulate a sequencing read or subpath, SCRaMbLE-SIM cuts the original chromosome at two positions. The first position is chosen randomly in the chromosome (uniform distribution), and the second position is chosen L LUs distant from the first cut, where L is the read length. 
The read length L is determined by a "Discretized Truncated Normal Distribution" (DTND). As a default, we used a mean read length of 8 LUs (~23.2 Kb) and a sigma of 3 LUs (~8.7 Kb).

All the essential functions for simulating long sequencing reads are in the script SCRaMbLE_DNA_simulation.py. This can be called using the function **DNA_extraction_coverage**.

```py
DNA_extraction_coverage(syn_chr, coverage: int, reads_len=8, sigma=3, circular=False)
```

I am also working on the script simulate_reads.py which will take as input a genome in fasta file and output simulated long reads in a fasta file.

### <a name="CN"></a>Genome composition and segment copy number analysis
Using simulated chromosomes, it is possible to investigate the fate of LUs during SCRaMbLE evolution. This means studying which segments are easier to be lost or to be conserved
The segment copy number after many SCRaMbLE simulations can be plotted with the function **SCRaMbLE_SIM_LU_CN** in the Mapping_coverage_MM.py script.

```py
SCRaMbLE_SIM_LU_CN(syn_chr, events=100, simulations=1000, essential=[], CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], max_CN=5)
```

By analysing these CN plots from simulated data, it might be possible to predict the results of future biological experiments, such as which parts of the chromosomes are more likely to be deleted or duplicated.
Additionally, studying the differences between simulated and biological data might give clues about the potential fitness effect of LUs. For example, LUs with a higher-than-expected CN might have a positive fitness effect, while low CN might indicate a negative fitness effect.

### <a name="population"></a>SCRaMbLE population SCRaMbLE population simulator
The SCRaMbLE population simulator can be used to study the SCRaMbLEd population structure.

So far, we have modelled and simulated evolution through a single genome/chromosome that keeps accumulating SCRaMbLE events.
However, a different approach consists of modelling and simulating the dynamics of an entire SCRaMbLE population where each cell has a different number of SCRaMbLE events (SEs) and cells share common ancestors and, therefore, some SEs.

Here we simulated a SCRaMbLE experiment where the Cre recombinase gene is placed under a daughter-specific promoter, and it is expressed only in the daughter cell but not in the mother cell (Lindstrom & Gottschling, 2009).
Therefore, at each generation, the cells replicate, and only the daughter cells undergo an x number of SCRaMbLE events.

```py
SCRaMbLE_population2(syn_chr, initial_cells=1, number_replication=8, events_for_replication=1, essential=[], mu=0, sigma=10, CEN=[], probability=[3, 2, 2, 1])
```

### <a name="mortality"></a>Mortality rate analysis of SCRaMbLE in a population
As the synthetic genome of a cell is SCRaMbLEd, there is a probability that an essential gene gets deleted and the cell dies, and we call this probability mortality rate. Indeed, this probability depends on the number of essential LUs, their position in the genome, the Cre activity, and the relative probability of deletions compared to other rearrangements.

To calculate the SCRaMbLE event mortality rate in chromosomes, you can use the function **percentage_of_chrs_with_essential_LUs_deleted**, which simulates 1000 times one SCRaMbLE event in the chromosome. Then it calculates the percentage of chromosomes with at least one essential LU deleted.

```py
percentage_of_chrs_with_essential_LUs_deleted(Chr, simulations=1000, Number_events=1, essential=[], mu=0, sigma=10, CEN=[], probability=[3, 2, 2, 1]):
```

Use the function **simulate_SCRaMbLE_pop_check_survival_rate** to simulate a SCRaMbLE population and calculate at every generation what is the mortality rate and how many cells are alive and dead.

```py
simulate_SCRaMbLE_pop_check_survival_rate(syn_chr, initial_cells=1, number_replication=8, events_for_replication=1, essential=[], mu=0, sigma=10, CEN=[], probability=[3, 2, 2, 1])
```

### <a name="chr_size"></a>Chromosome size evolution analysis
SCRaMbLE-SIM can be used to simulate many SCRaMbLEd chromosomes and perform some statistics about chromosome length evolution.

In the script SCRaMbLE_simulation_chr_len_events_easy.py there are all the functions created to perform this analysis. With the function **plot_SCRaMbLE_chr_len** we averaged the chromosome length during a 1,000 SCRaMbLE events evolution.

```py
plot_SCRaMbLE_chr_len(syn_chr, events=15, simulations=100, essential=[], CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], file_name="", SD=False)
```

However, we also study how the chromosome length during a SCRaMbLE evolution is influenced by several parameters like initial chromosome length, number of essential LUs, and relative probability of deletions and duplications.

```py
# This function plots the chromosome length over a SCRaMbLE evolution using different initial chromosome length.
chr_len_range_SCRaMbLE(events=15, simulations=100, essential=[], CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], file_name="")

# This function plots the chromosome length over a SCRaMbLE evolution using different number of essential LUs.
chr_len_essential_range_SCRaMbLE(syn_chr=50, events=15, simulations=100, CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], file_name="")

# This function plots the chromosome length over a SCRaMbLE evolution using different probabilities of deletions and duplications.
chr_len_probabilities_range_SCRaMbLE(syn_chr=50, events=15, simulations=100, essential=[], CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], file_name="", LOG=True)
```

## <a name="help"></a>Getting help
You can raise an issue on the issue page if you run into any problems or have any additional questions or requests.

## <a name="citation"></a>Citation
If you use this software, please consider citing our manuscript.

## <a name="limitations"></a>Limitations
SCRaMbLE-SIM was mostly tested on synthetic yeast genomes. For all other WT genomes which do not contain loxPsym sites, the recombination junctions can be chosen at random, or by the user to have equal size segments or chosen arbitrarily.

SCRaMbLE-SIM allows only genome evolution through genomic rearrangements and does not allows SNPs and INDELs.

At the moment, simulated reads have no errors inside and therefore might be different for this aspect to nanopore sequencing reads. Moreover, all simulated reads start and end at one recombination junction.

For the read length distribution, we chosen a "Discretized Truncated Normal Distribution" (DTND). However, this is just one of the three most common distributions for nanopore read length (Yu Li et al., 2018).
