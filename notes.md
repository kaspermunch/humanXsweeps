# TODO

- Copy script from argweaver-clues project and add to workflow_clues.py to summarize all clues info from files that finished.
- Make notebook that shows any peaks and shows allele frequency trajectories of any sweeps.


- Run Slim simulations in epoque
- Run notebook with dist to archaic
-

- [x] Map decode rec map to hg19
- [ ] NJ version of Tishkoff plots (check that placement of Ust ishim stays in ECH clade).
- [ ] Plot or report of paired Admixture proportions at each peak.
- [ ] See if any regions do not map to chrX in hg38.
- [ ] Change `prop_swept` to `min_prop_swept` in extended region data sets (do that on a clean up-to-date git version).
- [ ] Show peak plots in SOM using different clique sizes.
- [ ] Analyze properties of ECHs from truncated simulations.
- [ ] Run ARGweaver on each peak region using recombination map and Tennessee demography.
- [ ] Run CLUES on each region.

Address what the time since sweeps have done to the haplotype structure. How much is recombination expected to have broken down the sweep. How much more dramatic would they have needed to be to look like what they do now?


- Diversity:
    - [x] Exclude PARs from analysis (in analysis_globals) in both simons and g1000. I.e filter (start > 2699520) & (end < 154931044).
    - [x] Make Figure 1 and related figures from `male_dist_data_chrX` and not `all_male_dist_twice`
    - [x] Make sure I filter sweep_data using `not_missing` everywhere I should.
- Derived variants across X and autosomes:
    - [x] Counts of derived segregating variants
    - [x] Add Moises analysis to SOM.
    - [x] Counts of fixed (chimp and ancestral g1000)
- Calling of sweeps in Simons:
    - [x] Check strange haplotypes with no derived when using ancestral Tishkoff plot. (Missing data?)
    - [x] Look at and rerun Ust' Ishim code.
    - [x] Add dashed line with cutoff distance to the Tishkoff plots in papers.
    - [x] Test new way of calling sweeps and make sure it works better for `_af` calling.
    - [x] Run Tishkoff plots for all regions. Also make set with with missing data, and set showing admixture.
    - [ ] See if any regions do not map to chrX in hg38.
    - [ ] Change `prop_swept` to `min_prop_swept` in extended region data sets (do that on a clean up-to-date git version).
- Depletion of admixture in sweeps:
    - [x] Plot relative reduction in admixture against proportion swept. (Proportional reduction in admixture vs proportion swept)
    - [x] Plot distribution of admix_prop in ECHs.
    - [x] Make a new 100kb data frame with Laurits complete calls and merge with old one.
    - [x] Replicate analysis and make new cool plot.
- Distance to archaic humans (maybe):
    - [ ] Add Vindija Neanderthal to analysis. 
    - [ ] And the analysis to SOM.
- Analysis of 1000 genomes:
    - [ ] Get admix files from Laurits and make admix masking work for 1000 genomes.
    - [ ] Make a g1000_clique_calling.py
    - [ ] Exclude PARs from analysis as in Simons.
    - [ ] Get admixture data. 
    - [ ] Add admixture segment mask ("_af" files). 
    - [ ] Call sweeps on X and autosomes and compare them.
    - [ ] Do Africans have sweeps? suggesting low pi in regions is not BGS but also old sweeps. 
    - [ ] Do all African populations have lower diversity in regions?
    - [ ] Build 100kb data frame with admixture and replicate depletion analysis.
    - [ ] Make plot of admixture peaks on chr2.
    - [ ] Make sure admixture peaks are not artefacts (from over-collapsing of amplicons?), but maybe do not present them in this paper. See also note below.
    - [ ] See if admix depletion mirror population specific sweeps.
    - [ ] Argue that different heights of peaks could be ascribed to drift, but that qualitative differences shared by several populations cannot be the results of drift.
    - [ ] Do Fst analysis
    - [ ] Redo HapDAF run with physical distances
    - [ ] Do hapdaf analysis.
- Arguments for limits on neutral effects on X:
    - [ ] Is African diversity lower in All g1000 pops (in relation to start X/A for Pool-Nielsen)?
    - [ ] Find non-swept regions that match swept regions in terms of pi and recombination rate (that should also be peaks if peaks are driven by neutral processes)
    - [ ] Is a strong neutral effect compatible with no sweeps on autosomes?
- Arguments for how much admixture related BGS can reduce Ne:
    - [ ] Back-of-the-envelope argument for SOM
- Simulations:
    - [ ] Change to using cliques in simulations.
    - [ ] Set up a simulation truncated at 45,000 years and call ECHs using a pairwise distance cutoff of 2 * 4.2e-10 * 10000 = 8.4e-06).
    - [ ] Figure out for sure if maternal rate is per generation or per meiosis.
    - [ ] Do TMRCA half on slim simulations. 
    - [ ] Make sure simulated sweeps are only used if they complete (partially of fully). Make them conditional on completion.
    - [ ] Simulate sweeps with different coefficients.
    - [ ] Show ECH peaks for neutral and different types of sweeps.
    - [ ] 
- Estimation of selection coefficients:
    - [ ] ...
- Overlap with sites from Jarvis 2012:
    - [ ] Do Fisher exact test.
- Neutral simulations of allele frequencies:
    - [ ] 	Double-check simulations.
- Power and Ne for admixture inference:
    - [ ] Add Laurits simulations to SOM.
- Origin of population contributing swept haplotypes:
    - [ ] Add through explanation to SOM.




# Why strong peaks of admixture and why do some sweeps

If this is not an artefact high frequencies of admixture must be adaptive introgression - potentially as drive. One explanation for the overlap to sweeps is that archaic humans swept first and then separate humans swept later. But if such sweeps of Neandertal haplotypes are as recent as the sweeps we find, then the haplotypes should also be long. Also the sweeps should be strong given the small amount of time available. I should check if the contributions in each individual are from the same or different segments. An alternative scenario is that humans swept Neanderthals and then Neanderthals swept humans, carrying some flanking neanderthal along on the sweeping haplotypes. In this case the enters of the sweeps should be without admixture and the extra admixture should be in the flanks only.

Is it possible that Altai distance would not reveal exchange it if it is not the neanderthals we met?








# Ideas from EMBL conference
- Maybe try S-prime as backup for Laurits method to make absolutely certain that the regions are without admixture.
- Maybe look at non-ECH haplotyes at peaks to see if the 100kb distance to ECHs is small (suggesting that the peak is higher and more pointy, but that we just cannot see the sharp peak with our method). This would somewhat explain that swept regions have lower admixture than swept regions.
- Work on alternaive means to make sure the Ust Ishim shares haplotypes. Visualize In tishkoff plots. 
- Make it clear in the paper that our observation is admixture free ECHs in 10 Kyr.
- Make clear that we are not explaining all past observations of missing admixture on X.
- Compute X/A with and without regions and see if you can make it fit both X and Autosome distribution of distances or some other summary.
- Redo on 1000 genomes
- Look at LD in Africa to argue that the regions are subjected to drive too. Maybe also look at ARGweaver results among africans, e.g. rel_TMRCA_half or rel_TMRCA_75
- Simulate hybrid incompatibility, maybe as under-dominance in SLiM.
- It should be made clear how our observations relate to other reports of missing admixture on X. 
- Akey will send his segments from IBDmix
- The paper should not conclude the regions are sweeps until I have shown the admixture and sweep results.
- I need a metric to replace ECHs that is unambiguous and that only represents clades that form quickly. I.e. coalescences 50,000 and not just between present and 50,000.
	- Maybe run Relate to find single 100kb windows with
- Make sure it is not a problem that Laurits mask repeats and I do not, i.e that candidate regions are not special in terms of repeat content.
- Look at winiHS
- Run Relate and do sweep stat trick to estimate probability of  n lineages coalescing into n-10 in a time interval
- Show mean_dist for each 100kb window with peaks
- Look at regions that does not map to X in hg38




ECHs may be likely given some scenario... but that is not what carries our argument.
What carries it is that ECHs are without admixture.


I should actually compute X/A heterozygosity in African *outside* our regions to see if it is also 0.66 there. Assuming our regions take up p=0.1 of the chromosome.




The pair-wise distance plots do not show weather the excess of many small pair-wise distance is due to single large clades at each site, or one clade (as with a sweep).

Maybe ECH In simulations are more than One haplotype. 

Maybe try to call ECHs with a higher cutoff e.g. 50%

Maybe only do recombination mean for the 17 regions above 50%

Redo frequency simulation with the right parameters. 

What we want is to compute how likely it is that half of all lineages finds a common ancestor in less than 10,000 years given a harmonic mean bottle population size. I need to take into account that more than half may remain (using Asger’s trick) and average over the probabilities of different numbers of sample sizes at the onset of the bottle. I can use the harmonic mean pop size until the bottle. 

I can simulate that sweeps must have happened before 45kyr by splitting into three pops at that time. 


Pop structure. Split into five populations after bottleneck. 

Asger’s stats to compute probability of seeing 42 of lineages coalesce into one ancestor in 10000 years or in the time between OOA and population structure. 


Maybe 100kb windows in autosomes. 

Do simulation with mean recombination rate and no local Ne reduction for comparison to Figure 1

- While it is difficult to completely rule out a scenario in which neutral processes may give rise to such regions, we cannot explain why these regions are without admixture. 
- It is possible that one individual was not admixed and gave rise to a clade, but that is unlikely to have happened 17 independent times. 
- Strong negative selection would remove admixture but can not explain why admixture is only absent from high frequency haplotypes while it remains in the other individuals. It only works if promotion of individual haplotypes is linked to the absence of admixture. Negative selection would demote individual haplotypes with admixture (and secondarily promote all other haplotypes) but not promote *individual* haplotypes without it. 
- Demography: an X/A Ne ratio of 0.66 produces an X/A pi of 0.55, matching West-Eurasians that are (mostly) only affected by the OOA bottleneck.
- Use region/global pi in Africa (0.73) as proxy for how much Ne is lower in our regions. This reduces the X/A ratio to 0.66 * 0.73 = 0.48.

- Recombination: use chromosome wide recombination rate and I see no peaks. 
- A male bias of 2/1 (66% males) reduces the x/auto ratio from ~0.67 to 0.6.

>>> [0.67/0.75 * ((m * 1) + (1-m) * 2) / ( (m*2)+(1-m)*2) for m in [0.5, 0.66, 0.7, 0.75]]
	[0.67, 0.5985333333333334, 0.5806666666666668, 0.5583333333333333]

So I simulate with male bias of 2to1 and 4to1.




## Simulations to refute neutral theory:
- make sure sweeps have time to fix.
- compute actual X/auto heterozygosity in simons individuals, for use in simulations
Arguments for simulation parameters:
- Get phastcons in the regions.
- Show with pool Nielsen plot that a bottleneck is not expected to exacerbate variation in local Ne, making it the same as on autosomes: 
- Compare to empirical distance distribution.
- Do X and Auto neutral standard with a reduced scaled down demography.
	- This is to see when I start to see peaks without sweeps (It hopefully needs to unrealistically low).
	- Compare to empirical distance distribution to see (hopefully) that it is totally off.
	- See if a mixture of weak and super-strong bottleneck fits.
- Figure out why partial sweeps are above 50%.
- Refactor peak calling code and make `slim_trees.py` call peaks on each simulation.
- Plot 75% or 50% size against selection coefficient.
- How about TRMCA half? Is it still valid to compare swept and non-swept regions if they have different demographies?







# Links

https://www.simonsfoundation.org/simons-genome-diversity-project/
https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data/
https://docs.cancergenomicscloud.org/v1.0/docs/simons-genome-diversity-project-sgdp-dataset



# For later
- See if extra archaic admixture in Oceania swept regions is Neanderthal. (total length of segments with most snps from each archaic). The trough of Neanderthal should be more narrow than for denisova if it is secondary sweeps. 

# Old notes on genes

Among other notable genes overlapping sweeps are:… (needs to be updated after new sweep calling)
XIST/TSIX: Selection on skewed inactivation?
PIM2: (go:Spermatogenesis) kinase is one of the three highly conserved Pim family members which are known to be involved in cell survival and cell proliferation. (https://www.ncbi.nlm.nih.gov/pubmed?cmd=search&term=9804974&dopt=b) This gene encodes a protooncogene that acts as a serine/threonine protein kinase. Studies determined the encoded protein functions to prevent apoptosis and to promote cell survival. 
Notable in the light that that  copy number gain of VCX attributed to NOA by an underlying mechanism  that  induced  cell  apoptosis,  inhibited  cell  proliferation  and  disrupted  cell  cycle  progression,  thus  VCX may be a potential biomarker for the susceptibility to NOA
AR: (androgen receptor)
AKAP4: The encoded protein is localized to the sperm flagellum and may be involved in regulation of sperm motility (also circRNA)
FAM46D: Expressed during meiosis (pattern 10)
MID1: Midline (also ubiquitous circRNA)

Mental retardation / cognitive disability / psyciatric in OMIM

# Out of africa


# Precautions
- Rec rate in swept regions
- Consider removing B samples and these that Reich also removed: "5 samples based on missing X chromosome data in an initial processing, for themselves or a second sample: S_Finnish-1, S_Finnish-2, S_Mansi-1, S_Mansi-2, S_Palestinian-2"
- Mask out ampliconic genes before calling sweeps to make sure sweeps are not truncated by high pi from over-collapsing of ampliconic genes.
- Check afilt sweeps. 
- Use filtering level 9


# Hypotheses about ghost population
- What if there was a big standing variation of haplotypes in Africa. The first wave carried most of these but preferentially those from East Africa, if that was where this wave left Africa. The reason we see Africans as part of sweeps is then simply because these are the ones that share the particular swept haplotype. We know from Julien's study that the low ILS regions have lower recombination rate so it may not be unreasonable to think that these haplotypes could the shared across the set of African populations that are closet related to non-Africans. It then also makes sense that we do not see any individuals part of sweeps from the populations furthest from non-Africans (Ju_hoan_North, Mbuti, Baika, BantuTswana). We of curse need to check that the africans in sweep clades are  actual shared long haplotypes and not just random association expected from their closer relationship. If the hypothesis is true we would expect that most sweeps in non-africans share that haplotype with an African. Try to plot colored boxes where an African is part of a sweep and see if it covers most sweeps. This scenario would not require a stable high Ne low N first wave population, but simply that the first wave was bottlenecked and expanded to large N and that subsequent migrations (waves) together represented lots of variation. It could also explain why sweeps also remove Denisova admixture if sweeps are a *second* wave spreading across non-Africa. It also nicely explains why most sweeps are from a single or a few haplotypes (If sweeps came from a big standing variation in some ghost population, there should be many more haplotypes sweeping). Could be sweeps are not one single haplotype, but fixation of a mix of all the haplotypes in the admixture free population. 

Relationship with East Africa:
1. Kenyan variants swept after leaving Africa and was not swept in Africa.
    * Pi of Kenyan populations should not be low
    * male haplotypes of Kenyan populations should be identical to swept ones.
2. Sweeping variants also swept Kenyan populations after sweeping Middle East,
    * In this case the swept Kenyan populations should be more related to swept non-africans than other east African populations.
    * Maybe eastern Africa was part of one population genetic continuum with the middle east at the time of the sweeps.
    * Maybe the main wave was already bottlenecked in Kenya so that the second wave would also sweep Kenyan populations.
3. Variants were swept in East Africa and then later in Middle East.


Two alternatives:

1. Swept in the Middle East and spread from there along with rest of genome:
    - Supported by the fact that sweeps are partial and not shared by all males in each population (populations are a mix of populations that were swept and populations that were not).
    - Would require that sweeps were after Denisova admixture.
2. Rest of genome spread first and then the sweeps spread:
    - Supported by lack of Denisova admixture if sweeps postdate first wave.
    - Would make sense if the selection coefficient of drivers are frequency dependent: the sweeps do not complete if the selection coefficient goes down as the driver displaces the Y haplotype it drives against. 

Can they be distinguished?


The larger enrichment in Jegou et al. when controlling for heterozygosity.



# On meiotic drive
Problem: if the x is swept by segregation distorting haplotypes then the distortion must be strong, like >5%. But this would reduce the effective population size to about half in roughly ten generations. I.e. a compensatory mutation would need to arise in a very short time span. However, if pairs are monogamous then the excess females would just be lost in each generation and the ratio would stabilise at e.g. 55%

Potential explanation of how variants from populations with large standing variation can sweep through populations with small sanding variation. This can happen if transmission distortion is not controlled by the variant on X but which variants that are on the Y too. I.e. The transmission distortion of an X variant is determined by its interaction with variants on Y. With large standing variation any allele will only have a small frequent of partners with which it produces distortions. So averaged over the genetic backgrounds the selection coefficients are small. However, if a large populations meets a small population, then some/one of the many alleles will have a distortion advantage with one of the very few Y alleles in the small populations. Because the diversity of the "Y interaction background" is then small, the variation in background-averaged selection coefficients will be large.

# Devil's advocate:
- sweeps are just bottleneck
- How do I counter the argument that the high frequency of one divergence time is no just because that is where the bottle neck was? I guess the argument is that it should then apply randomly to the chromosome, and not to particular chromosomal regions. Also it should apply to the entire genome. Show that there are only very few similar depletions of diversity on the chr7 (do we expect the distribution of distances between pseudo haploids to be the same as that between true haploids? I guess they are only the same if sweeps are absent or completed. If sweeps are partial, then individuals heterozygous for the the sweep will show half the reduction in diversity in each pseudo haploid. Maybe it is most fair to compare pair wise distances for female X and chromosome 7  You would need to do all distances for x and 7 then.). 
- Bottleneck + Pool-Nielsen: If the coefficient of variation is expected to be the same on X and A (irrespective of global chromosome Ne), do we then expect pool-nielsen to create more variation across X? Also a low region would need to be one haplotype without recombination since the bottleneck. 
- BGS: That would apply to africa too then, and it does not.
- Heterosis: second wave had deleterious variants in higher freq. than first wave and admixed haplotypes rose in frequency due to heterosis. Not likely considering the time span for these sweeps.
- Interaction between sweep calling and admixture calling...
- If admixture is *after* sweeps we will not call sweeps in regions with admixture, and hence produce an apparent depletion of admixture in swept regions that is an artefact. In principle Neanderthal admixture could have happened after sweeps. If so, we would need to explain why some regions do not have less. That could be because of introgression barriers. If Denisova admixture is after sweeps then the same would need to apply. 
- If sweeps are *before* admixture, we can call sweeps just the same (some admixed segment would just be swept along).
- If sweeps are after neanderthal admixture before denisova admixture we should be able to call viewer sweeps in oceania, but that does not seem to be the case.

- Low admxiture is not characteristic of low pi regions. It is admixture that is characteristic of non-low pi regions. So it is possible that admixture introduces variation in non-pi regions thus inducing the observfed relation to pi. The effect would have to be stronger on X than on autosomes and is expected to be if admixture proportions arre the same on X and autosomes. The proportional incraese in pi in regions of admxiture versus regions without is: $$\frac{(1-\alpha)4N + \alpha(4N + 2T)}{4 N_e} = \frac{4N + 2 \alpha T}{4N} = \frac{2N + \alpha T}{2N}$$ If split is 500,000 years = 20,000 generations, $\alpha$ is 0.05, N is 10000 for autosomes and 5000 for chrX then the reduction is ~10% for chrX and ~5% for autosomes. However, $\alpha$ is an order of magnitude smaller ~0.005, implying that the relative increase would only be ~1%.
- The density of sites subject to BGS could be so high in the low pi regions that it would have resisted introgression no matter if they were located on X or autosomes. The logistic regression plots suggest that this is not the case.
- Admxiture with hybrid incompatibility could induce stronger BGS and cause depression of the low pi regions. However, this effect is likely not very strong as admix proportions are so low.

# For later
- The next paper can be a detailed account of how sweeps distribute across non-Africa.
- See if there is balancing selection in African windows that are swept in non-Africans.
- Maybe use argweaver to compute a statistic to measure balancing selection 
- Re-read argweaver paper and check that I handle discretization correctly.
- Explaining low X/A: Try to do X/A ratio vs Ne only on surely_not_swept regions.


## Running ARGweaver
Find a way to get discretization without hardcoding it here best way would be to pass a times file to the template as argument. I should make the argweaver template keep the times file and make it an output file. That way the --discretization argument to tree_statistics.py could be the times file and not the long list of values.
    
Tree statistics:

1. Sweep score
2. TMRCA (my_tmrca): time to TMRCA of tree.
3. TMRCA halftime (my_tmrca_half): time until at a clade consist of at least 50% of sample.
4. Coalescence half time (coal_half): time until at least half the sample is part of any subtree (have coalesced). Same as TMRCA_half but not requiring that they are in a single subtree.
5. Sweep node count (SwNC): Number of nodes with a short distance up and at least one sort distance down. Top node with no parent is not included.
6. Short branch enrichment (ShBE): Count how many of the shortest branches sort branches that together make up no more than e.g. 10% of total tree length.
7. Largest short branch tree (LShBT): Find the largest set of connected branches that are all smaller or equal to some maximum length.

Traverse tree if a node has only non-African descendants, count it and do not traverse further down below that node.
The number of African ancestors for non-africans is that count plus one (because the root counts as one too)

- Date sweep using tmrca_90


 
# argweaver ideas
 
1. Consider a clade root node.
2. Length of discretisation interval: t
3. Number of live lineages in the graph at beginning of interval: n
3. What is the expected time for coalescence of the lineages of the clade entering the interval (k).

it is  not just sum_{i=0}^{n-k} (n-i)(n-i-1)/2 because that does not take into account that the coalescences interval that are not part of the clade...
 
 
Maybe plot TMRCA along with pi to see that regions where some populations have low pi do not have lower than expected tmrca.

Make a GWF task that samples populations or regions separately. 
 
The clade that represents the time interval with the fastest rate of coalescence into one common ancestor (taking into account that this happens more often by chance if there are many life nodes)

Sum_i^{k-1} (k-i)(k-i-1)/2 

Expected depth of sweep clade / decided by depth of total tree from there on
 
 
Idea: Maybe draw a line between all pairs of populations that share a sweep
 
Maybe do an MSMC on each population to see if any bottlenecks fit our idea that this explains that sweeps are only outside africa.
 
### PCA of SNPs vs Pi
 
1. Do a PCA of individual SNPs to recover Nike pattern.
2. Do the same PCA wth all homozygotes encoded as 0 and see if Nike pattern is lost.
3. If it is lost then do PCA on Pi to on larger and larger bins to see when it is recovered.
 
If it turns out that there is some true signal in pi that is not just SNP information, then try to do PCA on pi for subsets with different regional pi. Try to run one for only mega base regions with high pi in all populations and see if the Nike pattern is lost (because there are only few sweeps)
    
 



# Gene info



**Entrez Gene Summary for AKAP4 Gene:** 
The A-kinase anchor proteins (AKAPs) are a group of structurally diverse proteins, which have the common function of binding to the regulatory subunit of protein kinase A (PKA) and confining the holoenzyme to discrete locations within the cell. This gene encodes a member of the AKAP family. The encoded protein is localized to the sperm flagellum and may be involved in the regulation of sperm motility. Alternative splicing of this gene results in two transcript variants encoding different isoforms. [provided by RefSeq, Jul 2008]

**GeneCards Summary for AKAP4 Gene:**
AKAP4 (A-Kinase Anchoring Protein 4) is a Protein Coding gene. Diseases associated with AKAP4 include Retinitis Pigmentosa 70. Among its related pathways are Activation of cAMP-Dependent PKA and Signal transduction_PKA signaling. Gene Ontology (GO) annotations related to this gene include protein kinase A binding. An important paralog of this gene is AKAP3.

----

**Entrez Gene Summary for NUDT10 Gene: **
This gene is a member of the nudix (nucleoside diphosphate linked moiety X)-type motif containing family. The encoded protein is a phosphohydrolase and may regulate the turnover of diphosphoinositol polyphosphates. The turnover of these high-energy diphosphoinositol polyphosphates represents a molecular switching activity with important regulatory consequences. Molecular switching by diphosphoinositol polyphosphates may contribute to the regulation of intracellular trafficking. In some populations putative prostate cancer susceptibility alleles have been identified for this gene. Alternatively spliced transcript variants, which differ only in the 5' UTR, have been found for this gene. [provided by RefSeq, Feb 2015]

**GeneCards Summary for NUDT10 Gene:**
NUDT10 (Nudix Hydrolase 10) is a Protein Coding gene. Diseases associated with NUDT10 include Autoimmune Disease Of Endocrine System. Among its related pathways are Metabolism and Inositol phosphate metabolism (REACTOME). Gene Ontology (GO) annotations related to this gene include hydrolase activity and inositol diphosphate tetrakisphosphate diphosphatase activity. An important paralog of this gene is NUDT11.

----

**GeneCards Summary for CXorf67 Gene:**
CXorf67 (Chromosome X Open Reading Frame 67) is a Protein Coding gene. Diseases associated with CXorf67 include Endometrial Stromal Sarcoma and Endometrial Stromal Nodule.

----

**Entrez Gene Summary for WNK3 Gene:**
This gene encodes a protein belonging to the 'with no lysine' family of serine-threonine protein kinases. These family members lack the catalytic lysine in subdomain II, and instead have a conserved lysine in subdomain I. This family member functions as a positive regulator of the transcellular Ca2+ transport pathway, and it plays a role in the increase of cell survival in a caspase-3-dependent pathway. Alternative splicing results in multiple transcript variants. [provided by RefSeq, May 2010]

**GeneCards Summary for WNK3 Gene:**
WNK3 (WNK Lysine Deficient Protein Kinase 3) is a Protein Coding gene. Diseases associated with WNK3 include Syndromic X-Linked Intellectual Disability Siderius Type and Lung Large Cell Carcinoma. Among its related pathways are Diuretics Pathway, Pharmacodynamics and Ion channel transport. Gene Ontology (GO) annotations related to this gene include transferase activity, transferring phosphorus-containing groups and protein tyrosine kinase activity. An important paralog of this gene is WNK2.

**UniProtKB/Swiss-Prot for WNK3 Gene:**
WNK3_HUMAN,Q9BYP7
Serine/threonine kinase which plays an important role in the regulation of electrolyte homeostasis, cell signaling, survival and proliferation. Acts as an activator and inhibitor of sodium-coupled chloride cotransporters and potassium-coupled chloride cotransporters respectively (PubMed:16275913, PubMed:16275911, PubMed:16357011). Phosphorylates WNK4. Regulates the phosphorylation of SLC12A1 and SLC12A2. Increases Ca(2+) influx mediated by TRPV5 and TRPV6 by enhancing their membrane expression level via a kinase-dependent pathway (PubMed:18768590). Inhibits the activity of KCNJ1 by decreasing its expression at the cell membrane in a non-catalytic manner.

WNK3 is associated with the GO term: GO:0043066	**negative regulation of apoptotic process**

----

**GeneCards Summary for TCEAL6 Gene:**
TCEAL6 (Transcription Elongation Factor A Like 6) is a Protein Coding gene. An important paralog of this gene is TCEAL3.

**UniProtKB/Swiss-Prot for TCEAL6 Gene:**
May be involved in transcriptional regulation.

----

**Entrez Gene Summary for IL13RA2 Gene:**
The protein encoded by this gene is closely related to Il13RA1, a subuint of the interleukin 13 receptor complex. This protein binds IL13 with high affinity, but lacks cytoplasmic domain, and does not appear to function as a signal mediator. It is reported to play a role in the internalization of IL13. [provided by RefSeq, Jul 2008]

**GeneCards Summary for IL13RA2 Gene:**
IL13RA2 (Interleukin 13 Receptor Subunit Alpha 2) is a Protein Coding gene. Diseases associated with IL13RA2 include Malignant Glioma and Glioblastoma Multiforme. Among its related pathways are Cytokine Signaling in Immune system and Akt Signaling. Gene Ontology (GO) annotations related to this gene include signal transducer activity and cytokine receptor activity. An important paralog of this gene is IL5RA.

----

**Entrez Gene Summary for ACTRT1 Gene:**
This gene encodes a protein related to the cytoskeletal protein beta-actin. This protein is a major component of the calyx in the perinuclear theca of mammalian sperm heads, and it therefore likely functions in spermatid formation. This gene is intronless and is similar to a related gene located on chromosome 1. A related pseudogene has also been identified approximately 75 kb downstream of this gene on chromosome X. [provided by RefSeq, May 2010]

**GeneCards Summary for ACTRT1 Gene:**
ACTRT1 (Actin Related Protein T1) is a Protein Coding gene. An important paralog of this gene is ACTRT2.

----

**Entrez Gene Summary for H2AFB1 Gene:**
Histones are basic nuclear proteins that are responsible for the nucleosome structure of the chromosomal fiber in eukaryotes. Nucleosomes consist of approximately 146 bp of DNA wrapped around a histone octamer composed of pairs of each of the four core histones (H2A, H2B, H3, and H4). The chromatin fiber is further compacted through the interaction of a linker histone, H1, with the DNA between the nucleosomes to form higher order chromatin structures. This gene encodes a replication-independent histone that is a member of the histone H2A family. **This gene is part of a region that is repeated three times on chromosome X, once in intron 22 of the F8 gene and twice closer to the Xq telomere. This record represents the most centromeric copy which is in intron 22 of the F8 gene.** [provided by RefSeq, Oct 2015]

- H2AFB1 at chrX:154689080-154689596 - (NM_001017990) histone H2A-Bbd type 1
- H2AFB1 at chrX:154610428-154610944 - (NM_001017990) histone H2A-Bbd type 1
- H2AFB1 at chrX:154113317-154113833 - (NM_001017990) histone H2A-Bbd type 1

**GeneCards Summary for H2AFB1 Gene:**
H2AFB1 (H2A Histone Family Member B1) is a Protein Coding gene. Among its related pathways are Mitotic Prophase and Meiosis. Gene Ontology (GO) annotations related to this gene include protein heterodimerization activity. An important paralog of this gene is H2AFB2.

**UniProtKB/Swiss-Prot for H2AFB1 Gene:**
***Atypical histone H2A*** which can replace conventional H2A in some nucleosomes and is associated with active transcription and mRNA processing (PubMed:22795134). Nucleosomes wrap and compact DNA into chromatin, limiting DNA accessibility to the cellular machineries which require DNA as a template. Histones thereby play a central role in transcription regulation, DNA repair, DNA replication and chromosomal stability (PubMed:15257289, PubMed:16287874, PubMed:16957777, PubMed:17591702, PubMed:17726088, PubMed:18329190, PubMed:22795134). Nucleosomes containing this histone are less rigid and organize less DNA than canonical nucleosomes in vivo (PubMed:15257289, PubMed:16957777, PubMed:17591702, PubMed:24336483). They are enriched in actively transcribed genes and associate with the elongating form of RNA polymerase (PubMed:17591702, PubMed:24753410). They associate with spliceosome components and are required for mRNA splicing (PubMed:22795134).


One overlapping protein coding gene is:

**MAGEH1**: This gene belongs to the non-CT (non cancer/testis) subgroup of the melanoma-associated antigen (MAGE) superfamily. The encoded protein is likely associated with apoptosis, cell cycle arrest, growth inhibition or cell differentiation. The protein may be involved in the atRA (all-trans retinoic acid) signaling through the STAT1-alpha (signal transducer and activator of transcription 1-alpha) pathway. [provided by RefSeq, Aug 2013]


Another is **AR**:

**Entrez Gene Summary for AR Gene:**
The androgen receptor gene is more than 90 kb long and codes for a protein that has 3 major functional domains: the N-terminal domain, DNA-binding domain, and androgen-binding domain. The protein functions as a steroid-hormone activated transcription factor. Upon binding the hormone ligand, the receptor dissociates from accessory proteins, translocates into the nucleus, dimerizes, and then stimulates transcription of androgen responsive genes. This gene contains 2 polymorphic trinucleotide repeat segments that encode polyglutamine and polyglycine tracts in the N-terminal transactivation domain of its protein. Expansion of the polyglutamine tract from the normal 9-34 repeats to the pathogenic 38-62 repeats causes spinal bulbar muscular atrophy (SBMA, also known as Kennedy's disease). Mutations in this gene are also associated with complete androgen insensitivity (CAIS). Alternative splicing results in multiple transcript variants encoding different isoforms. [provided by RefSeq, Jan 2017]


Also associated with: GO:2001237, **negative regulation of extrinsic apoptotic signaling pathway**



X chromosome inactivation: lincRNA. X-inactivation center: XIST, TSIX, JPX, FTX.

Y RNA: Two functions have been described for Y RNAs in the literature: As a repressor of Ro60, and as an initiation factor for DNA replication. Mutant human Y RNAs lacking the conserved binding site for Ro60 protein still support DNA replication,[3] indicating that binding to Ro protein and promoting DNA replication are two separable functions of Y RNAs. Although Y RNA-derived small RNAs are similar in size to microRNAs, it has been shown that these Y RNA fragments are not involved in the microRNA pathway.[8]

Y RNAs are overexpressed in some human tumours and are required for cell proliferation[11] and small, microRNA-sized breakdown products may be involved in autoimmunity and other pathological conditions.[12] Recent work has demonstrated that Y RNAs are modified at their 3' end by the non-canonical poly(A) polymerase PAPD5, and the short oligo(A) tail added by PAPD5 is a marker for 3' end processing by the ribonuclease PARN/EXOSC10 or for degradation by the exonuclease DIS3L.[13] Since PARN deficiency causes a severe form of the bone marrow disease Dyskeratosis Congenita as well as pulmonary fibrosis,[14][15] it is possible that defects in Y RNA processing contribute to the severe pathology observed in these patients

From wikipedia







Entrez Gene Summary for PAGE4 Gene: This gene is a member of the GAGE family. The GAGE genes are expressed in a variety of tumors and in some fetal and reproductive tissues. This gene is strongly expressed in prostate and prostate cancer. It is also expressed in other male and female reproductive tissues including testis, fallopian tube, uterus, and placenta, as well as in testicular cancer and uterine cancer. The protein encoded by this gene shares sequence similarity with other GAGE/PAGE proteins, and also belongs to a family of CT (cancer-testis) antigens. The protein may play a role in benign and malignant prostate diseases. A related pseudogene is located on chromosome 7. Alternate splicing results in multiple transcript variants. [provided by RefSeq, Jan 2016]

GeneCards Summary for PAGE4 Gene: PAGE4 (PAGE Family Member 4) is a Protein Coding gene. Diseases associated with PAGE4 include Testicular Cancer and Prostate Cancer.







**Entrez Gene Summary for PAGE4 Gene:**
This gene is a member of the GAGE family. The GAGE genes are expressed in a variety of tumors and in some fetal and reproductive tissues. This gene is strongly expressed in prostate and prostate cancer. It is also expressed in other male and female reproductive tissues including testis, fallopian tube, uterus, and placenta, as well as in testicular cancer and uterine cancer. The protein encoded by this gene shares sequence similarity with other GAGE/PAGE proteins, and also belongs to a family of CT (cancer-testis) antigens. The protein may play a role in benign and malignant prostate diseases. A related pseudogene is located on chromosome 7. Alternate splicing results in multiple transcript variants. [provided by RefSeq, Jan 2016]

**GeneCards Summary for PAGE4 Gene:**
PAGE4 (PAGE Family Member 4) is a Protein Coding gene. Diseases associated with PAGE4 include Testicular Cancer and Prostate Cancer.







| GO term | Description | P-value |	FDR q-value | Enrichment (N, B, n, b) |
|:----|:----|:----|:----|:----|
| GO:0051054 | positive regulation of DNA metabolic process | 2.41E-4 | 1E0 | 6.55 (734,7,80,5) |
    

- DKC1: dyskeratosis congenita 1, dyskerin
- ATRX: alpha thalassemia/mental retardation syndrome x-linked
- FANCB: fanconi anemia, complementation group b
- PAK3: p21 protein (cdc42/rac)-activated kinase 3
- BRCC3: brca1/brca2-containing complex, subunit 3


BRCC3 sits in the same sweep region as the three copies of the H2AFB1 Gene that encodes an atypical histone H2A.

**Entrez Gene Summary for BRCC3 Gene:**
This gene encodes a subunit of the BRCA1-BRCA2-containing complex (BRCC), which is an E3 ubiquitin ligase. This complex plays a role in the DNA damage response, where it is responsible for the stable accumulation of BRCA1 at DNA break sites. The component encoded by this gene can specifically cleave Lys 63-linked polyubiquitin chains, and it regulates the abundance of these polyubiquitin chains in chromatin. The loss of this gene results in abnormal angiogenesis and is associated with syndromic moyamoya, a cerebrovascular angiopathy. Alternative splicing results in multiple transcript variants. A related pseudogene has been identified on chromosome 5. [provided by RefSeq, Jun 2011]

**GeneCards Summary for BRCC3 Gene:**
BRCC3 (BRCA1/BRCA2-Containing Complex Subunit 3) is a Protein Coding gene. Diseases associated with BRCC3 include Moyamoya Disease 4 With Short Stature, Hypergonadotropic Hypogonadism, And Facial Dysmorphism and T-Cell Prolymphocytic Leukemia. Among its related pathways are Metabolism of proteins and DNA Double-Strand Break Repair. Gene Ontology (GO) annotations related to this gene include metallopeptidase activity and obsolete ubiquitin thiolesterase activity.

**UniProtKB/Swiss-Prot for BRCC3 Gene:**
Metalloprotease that specifically cleaves Lys-63-linked polyubiquitin chains (PubMed:19214193, PubMed:20656690, PubMed:24075985, PubMed:26344097). Does not have activity toward Lys-48-linked polyubiquitin chains. Component of the BRCA1-A complex, a complex that specifically recognizes Lys-63-linked ubiquitinated histones H2A and H2AX at DNA lesions sites, leading to target the BRCA1-BARD1 heterodimer to sites of DNA damage at double-strand breaks (DSBs). In the BRCA1-A complex, it specifically removes Lys-63-linked ubiquitin on histones H2A and H2AX, antagonizing the RNF8-dependent ubiquitination at double-strand breaks (DSBs) (PubMed:20656690). Catalytic subunit of the BRISC complex, a multiprotein complex that specifically cleaves Lys-63-linked ubiquitin in various substrates (PubMed:20656690, PubMed:24075985, PubMed:26344097, PubMed:26195665). Mediates the specific Lys-63-specific deubiquitination associated with the COP9 signalosome complex (CSN), via the interaction of the BRISC complex with the CSN complex (PubMed:19214193). The BRISC complex is required for normal mitotic spindle assembly and microtubule attachment to kinetochores via its role in deubiquitinating NUMA1 (PubMed:26195665). Plays a role in interferon signaling via its role in the deubiquitination of the interferon receptor IFNAR1; deubiquitination increases IFNAR1 activity by enhancing its stability and cell surface expression (PubMed:24075985, PubMed:26344097). Down-regulates the response to bacterial lipopolysaccharide (LPS) via its role in IFNAR1 deubiquitination (PubMed:24075985).










