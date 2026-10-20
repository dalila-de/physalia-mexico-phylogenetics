# Prior genome-wide display assignments

These files document the existing E/W and vagrant labels retained at the author's
explicit request for the NCBI-based rerun. They are external display metadata,
not groups estimated from these marker trees or from GenBank annotations.

The recorded source is `Population_Paper_Phase1/admixture_likelihood_v1/results`.
The whole-assembly NGSadmix K2 component assignment was oriented toward the
eastern samples, with E assigned above 0.5 and W otherwise, checked against the
scaffold-10-excluded comparison. Prior vagrant flags were retained. These values
are not calibrated historical ancestry proportions. The source report and frozen
ancestry tables are retained here without modification.

`plot_provenance.json` describes the earlier display preparation and its original
artifact hashes. It is a historical record, not a checksum manifest for the new
trees. Current display metadata are in `analyses/inputs/genetic_annotations.tsv`,
with a fresh hash in the NCBI rerun provenance. The five Progreso specimens and
CH887 remain unassigned for E/W. No genome-wide inference was rerun here.

Church et al. (2025), [Population genomics of a sailing siphonophore reveals
genetic structure in the open ocean](https://doi.org/10.1016/j.cub.2025.05.066),
is the published population-genomic context; the specific retained assignments
are traced to the frozen source files above.
