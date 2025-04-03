import argparse

parser = argparse.ArgumentParser(description='Search for iterons around CRESS stem-loops in GFF files')

parser.add_argument(
    "--outputDir",
    help="Directory for output files (default: current directory)",
    default=".",
    type=str
)

parser.add_argument(
    "--inputFasta",
    help="Path to input FASTA file with all sequences",
    default="examples/test.fasta"
)

parser.add_argument(
    "--inputGFF",
    help="Path to associated input GFF file",
    default="examples/out.gff"
)

parser.add_argument(
    "--outputGFF",
    help="Path for output GFF file",
    default="examples/finaloutput.gff"
)

parser.add_argument(
    "--outputAnnotations",
    help="Identifiers to selectively preserve annotations",
    default="2 CRUISE"
)

parser.add_argument(
    "--minLength",
    help="Minimum iteron length",
    type=int,
    default=5
)

parser.add_argument(
    "--maxLength",
    help="Maximum iteron length",
    type=int,
    default=12
)

parser.add_argument(
    "--range",
    help="Number of base pairs around nona to search",
    type=int,
    default=65
)

parser.add_argument(
    "--rank",
    help="Use ranking system",
    type=bool,
    default=True
)

parser.add_argument(
    "--numberTopIterons",
    help="The number of iterons returned in rank order",
    type=int,
    default=5
)

parser.add_argument(
    "--maxScore",
    help="Maximum score allowed for iterons if rank = False",
    type=int,
    default=40
)

parser.add_argument(
    "--wiggle",
    help="Max difference between iteron length and distance allowed before score detriment",
    type=int,
    default=5
)

parser.add_argument(
    "--goodLength",
    help="The highest favorable iteron length",
    type=int,
    default=11
)

parser.add_argument(
    "--doStemLoop",
    help="Whether to annotate stem-loop repeats",
    type=bool,
    default=True
)

parser.add_argument(
    "--doKnownIterons",
    help="Whether to annotate known iterons",
    type=bool,
    default=True
)

parser.add_argument(
    "--maxDist",
    help="Maximum allowed distance between iterons",
    type=int,
    default=20
)

parser.add_argument(
    "--bestDist",
    help="Optimal maximum distance between iterons",
    type=int,
    default=10
)

parser.add_argument(
    "--scoreRange",
    help="Score range between outputted candidates (recommended 50)",
    type=int,
    default=50
)
