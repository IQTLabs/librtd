import librtd
import cli
import tables
import nimpy

func same_kmer_return_times*(x: string, k: Positive): Table[string, seq[int]] {.exportpy.} =
  librtd.sameKmerReturnTimes(x, k)

func pairwise_kmer_return_times*(x: string, k: Positive): Table[string, seq[int]] {.exportpy.} =
  librtd.pairwiseKmerReturnTimes(x, k)

func reverse_complement_return_times*(x: string, k: Positive): Table[string, seq[int]] {.exportpy.} =
  librtd.reverseComplementReturnTimes(x, k)

func return_time_distribution*(x: string, k: Positive, pairwise: bool = false, reverse_complement: bool = false): Table[string, float] {.exportpy.} =
  librtd.returnTimeDistribution(x, k, pairwise=pairwise, reverseComplement=reverse_complement)

proc main*(k: int, input: string, output: string = "stdout", reverseComplement: bool = false, pairwise: bool = false) {.exportpy.} =
  cli.main(k, input, output, reverseComplement, pairwise) 