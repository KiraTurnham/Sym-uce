# LSU reference panel

[`symbiodiniaceae_lsu.fasta`](symbiodiniaceae_lsu.fasta) contains 77 LSU
records intended for read enrichment and assembly, not species-level
identification:

- 75 sequences from the `[Family Symbiodiniaceae]` block in the LSU alignment
  `CurrentBiology_LSU_dinoflagellate_Orders.nexus`.
- Two GenBank 28S records: *Durusdinium trenchii* `KF740689.1` and
  *D. glynii* `KY131785.1`.

The alignment was downloaded by the user from Dryad as part of the dataset
associated with LaJeunesse et al. (2018), “Systematic Revision of
Symbiodiniaceae Highlights the Antiquity and Diversity of Coral Endosymbionts”
([article DOI](https://doi.org/10.1016/j.cub.2018.07.008);
[dataset DOI](https://doi.org/10.5061/dryad.1717129), CC0). Only entries in the
Symbiodiniaceae family block were included. Alignment gap characters were
removed; `?` characters were replaced with `N`; nucleotide ambiguity symbols
were retained. The NEXUS alignment labels are preserved in FASTA headers with
the genus added as a prefix.

The panel represents seven genera: *Cladocopium*, *Breviolum*, *Fugacium*,
*Durusdinium*, *Gerakladium*, *Effrenium*, and *Symbiodinium*. It is a broad
LSU enrichment reference, not a comprehensive or balanced census. In
particular, the C21/C40 labels occur in the separate COB alignment, which is
not used here. An unsampled *Cladocopium* species can still map to the
*Cladocopium* LSU records, but this panel cannot establish that it is C40 or
any other species.

See [`lsu_reference_sources.tsv`](lsu_reference_sources.tsv) for the source
label and citation associated with each FASTA record.
