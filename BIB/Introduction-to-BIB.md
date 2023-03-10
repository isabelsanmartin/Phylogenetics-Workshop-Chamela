# Introduction to the Bayesian Island Biogeography (BIB) Model in RevBayes

Instructor: **Isabel Sanmartin**

For a good introduction to the BIB model, you can read Sanmartín et al. (2008) or this book chapter [Sanmartin, 2022 ISME](https://books.google.com.mx/books?hl=en&lr=&id=7b5WEAAAQBAJ&oi=fnd&pg=PA27&dq=info:LWXmCk9U9N4J:scholar.google.com&ots=Sk8feF1yvv&sig=HP7U3SJzoSfGRyLrOXJk2DJZF4Y&redir_esc=y#v=onepage&q&f=false)

### Bayesian Island Model

The BIB model uses MCMC Bayesian Inference to estimate ancestral ranges and rates of biogeographic parameters alongside phylogenetic parameters, such as the tree topology, parameters of molecular evolution, and branch lengths; the input data are DNA sequences and tip distributions of the study species (Sanmartín et al. 2008).


Though they both use CTMC processes to model range evolution, BIB and DEC models are slightly different (Figure 1). BIB (Fig. 1a) implements a simpler character evolutionary model, in which ancestors can only occupy single areas (A or B) and range evolution along the branches if governed by a CTMC process with only one type of parameter equivalent to range switching or instantaneous dispersal; the Q matrix describes the instantaneous transition from one area as a jump dispersal event (p = A to B). At the speciation events in the phylogeny, the single-area ancestral range is inherited entirely and identically by the two descendants, in other words, there is no need to include a cladogenetic component in the BIB model because the ancestral range is not altered through speciation (Fig. 5a). 

![Figure0](figures/Figure0.png "Figure 0")*Parametric CTMC models in biogeographic inference. a. Standard CTMC model (BIB). b. Dispersal-Extinction Cladogenesis model (DEC)*


There are advantages and some disadvantages when comparing the BIB against the DEC model: hough they both use CTMC processes to model range evolution, BIB and DEC models are slightly different (Figure 5). BIB (Fig. 5a) implements a simpler character evolutionary model, in which ancestors can only occupy single areas (A or B) and range evolution along the branches if governed by a CTMC process with only one type of parameter equivalent to range switching or instantaneous dispersal; the Q matrix describes the instantaneous transition from one area as a jump dispersal event (p = A to B). At the speciation events in the phylogeny, the single-area ancestral range is inherited entirely and identically by the two descendants, in other words, there is no need to include a cladogenetic component in the BIB model because the ancestral range is not altered through speciation (Fig. 5a). 

<br>

Modeling dispersal as an instantaneous process without going through a widespread state may seem unrealistic but allows the BIB model to "borrow" the sophisticated machinery and statistical algorithms used in molecular models of nucleotide substitution; in fact, initial implementations of BIB used software routinely employed in molecular phylogenetics (Sanmartín et al. 2008; Lemey et al. 2009). 

<br>

In standard molecular models, nucleotide substitutions within a species DNA sequence are considered as instantaneous In-between demographic-level processes, involving increased allele polymorphism within gene trees, competition among mutations in terms of fitness, and rates of fixation differing between alleles (De Maio et al. 2015), are typically ignored.
Similarly, in the BIB model, the species is assumed to instantaneously change area relative to its current range, ignoring the intermediary population-level processes, such as changes in effective population size with migration, introgression, etc.
The BIB model can thus be appropriate to model scenarios in which areas are discrete entities isolated by dispersal barriers, so that migration to a new area effectively leads to speciation, in other words, the ancestor is not expected to maintain the widespread distribution for long, as in the case of founder effects in oceanic islands isolated by geographic barriers (Sanmartin et al. 2008), or in continental islands isolated by ecological barriers (Sanmartín et al. 2010).

<br>

However, the assumption of single-state ancestral ranges means that BIB is most useful to explore and test general patterns of geographic movement or dispersal; if the interest lies on inferring speciation modes or possible ways in which ancestral ranges are divided, BIB is not well suited.
Notice that constraining ancestors to single areas in the Q matrix does not imply that phylogenies with extant widespread species cannot be analyzed with BIB. As in molecular evolutionary models, these widespread terminals will be treated as sources of "ambiguity" in the BIB analysis: 50% of the time the MCMC chain will sample from one of the discrete states, and 50% from the other.
