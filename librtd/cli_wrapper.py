"""Return time distribution (RTD) calculation.

Takes input FASTA files and outputs a line-delimited JSON (.jsonl) file containing the RTD for each k-mer.
If no output file is specified, it will be written to stdout.
All log messages are written to stderr.

Usage:
rtd <k> <input> [<output>] [--reverse-complement|--pairwise]
rtd (-h | --help)
rtd --version

Options:
-r, --reverse-complement  Whether to compute distances to reverse complement k-mers
-p, --pairwise            Whether to compute the distances between every pair of k-mers
-h, --help                Show this screen.
--version                 Show version.
"""
from docopt import docopt
import nimporter
from .cli import main


def cli_wrapper():
    arguments = docopt(__doc__, version="Naval Fate 2.0")
    output = arguments["<output>"] if arguments["<output>"] is not None else "stdout"
    main(
        int(arguments["<k>"]),
        arguments["<input>"],
        output=output,
        reverseComplement=arguments["--reverse-complement"],
        pairwise=arguments["--pairwise"],
    )


if __name__ == "__main__":
    cli_wrapper()
