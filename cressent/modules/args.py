# arguments for run_cruise

import argparse

parser = argparse.ArgumentParser(description='Search for iterons around CRESS stem-loops in GFF files')

parser.add_argument(
    "--output",
    help="Path to the output directory (Default: working directory)",
    default=".",
    type=str
)

parser.add_argument(
    "--input_fasta",
    help="Path to input FASTA file with all sequences"
)

parser.add_argument(
    "--inputGFF",
    help="Path to associated input GFF file"
)

parser.add_argument(
    "--outputGFF",
    help="Path for output GFF file (default: finaloutput.gff)",
    default="finaloutput.gff"
)

parser.add_argument(
    "--outputAnnotations",
    help="Identifiers to selectively preserve annotations (default: 2 CRUISE)",
    default="2 CRUISE"
)

parser.add_argument(
    "--minLength",
    help="Minimum iteron length (default = 5)",
    type=int,
    default=5
)

parser.add_argument(
    "--maxLength",
    help="Maximum iteron length (default = 12)",
    type=int,
    default=12
)

parser.add_argument(
    "--range",
    help="Number of base pairs around nona to search (default = 65)",
    type=int,
    default=65
)

parser.add_argument(
    "--rank",
    help="Use ranking system (default = True)",
    type=bool,
    default=True
)

parser.add_argument(
    "--numberTopIterons",
    help="The number of iterons returned in rank order (default = 5)",
    type=int,
    default=5
)

parser.add_argument(
    "--maxScore",
    help="Maximum score allowed for iterons if rank = False (default = 40)",
    type=int,
    default=40
)

parser.add_argument(
    "--wiggle",
    help="Max difference between iteron length and distance allowed before score detriment (default = 5)",
    type=int,
    default=5
)

parser.add_argument(
    "--goodLength",
    help="The highest favorable iteron length (default = 11)",
    type=int,
    default=11
)

parser.add_argument(
    "--doStemLoop",
    help="Whether to annotate stem-loop repeats (default = True)",
    type=bool,
    default=True
)

parser.add_argument(
    "--doKnownIterons",
    help="Whether to annotate known iterons (default = True)",
    type=bool,
    default=True
)

parser.add_argument(
    "--maxDist",
    help="Maximum allowed distance between iterons (default = 20)",
    type=int,
    default=20
)

parser.add_argument(
    "--bestDist",
    help="Optimal maximum distance between iterons (default = 10)",
    type=int,
    default=10
)

parser.add_argument(
    "--scoreRange",
    help="Score range between outputted candidates (default: 50)",
    type=int,
    default=50
)
