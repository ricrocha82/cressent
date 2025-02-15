# ssDNA tool
a modular tool to help researchers to automatically annotate ssDNA contigs

# DB with all protein sequences
[Families](/fs/project/PAS1117/ricardo/ssDNA_tool/build_db/names_to_query_ssdna.csv) were collected from the ICTV website

Only the Pleolipoviridae Family (contains both dsDNA and ssDNA) were filtered by Genus names:
- Alphapleolipovirus finnoniense
- Alphapleolipovirus huluense
- Alphapleolipovirus samutsakhonense
- Alphapleolipovirus thailandense

Sequences that were not found in Refseq were directly download (if available) using the [Acession number](/fs/project/PAS1117/ricardo/ssDNA_tool/build_db/accession_not_found_refseq.csv) extracted from the reports found on [ictv report website](https://ictv.global/report):
- Adamaviridae
- Anicreviridae
- Draupnirviridae
- Endolinaviridae (not found)
- Gandrviridae
- Geplanaviridae (not found)
- Kanorauviridae (not found)
- Kirkoviridae
- Mahapunaviridae (not found)
- Ouroboviridae
- Pecoviridae

All sequences are here: "/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/families_aa"

## processing DB (for global comparison)
- the steps below were conducted by each family.
- clustered at 95% amino acid identity > 90% (`CD-HIT`)
- compare sequences using all-vs-all BLASTp with e-value 1e-5 (`diamond`)
- clustered using Markov cluster `MCL` algorithm - inflation rate at `1.5`
- clusters with annotated Rep and Capsid/VP will be concatenated and aligned using `MAFFT` with the `auto` parameter.
- trimmed using `Trimal` with gap threshold `0.15`
- to generate a global DB (all families) the concatenated clusters of all families were concatenated, aligned and trimmed using the same steps as above.

- families without Rep proteins: Bidnaviridae, Finnlakeviridae, Metaxyviridae, and Spiraviridae
- families without Cap proteins: Alphasatellitidae, Spiraviridae, and Tolecusatellitidae
