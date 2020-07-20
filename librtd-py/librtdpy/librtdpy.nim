import librtd
import tables
import nimpy

func kmer_indices*(x: string, k: Positive): Table[string, seq[int]] {.exportpy.} =
  librtd.kmerIndices(x, k)

func same_kmer_return_times*(x: string, k: Positive): Table[string, seq[int]] {.exportpy.} =
  librtd.sameKmerReturnTimes(x, k)

func dist_to_next_greater_index*(indicies1: seq[int], indices2: seq[int]): seq[int] {.exportpy.} =
  librtd.distToNextGreaterIndex(indicies1, indices2)

func pairwise_kmer_return_times*(x: string, k: Positive): Table[string, seq[int]] {.exportpy.} =
  librtd.pairwiseKmerReturnTimes(x, k)

func reverse_complement_return_times*(x: string, k: Positive): Table[string, seq[int]] {.exportpy.} =
  librtd.reverseComplementReturnTimes(x, k)

func return_time_distribution*(x: string, k: Positive, pairwise: bool = false, reverse_complement: bool = false): Table[string, float] {.exportpy.} =
  librtd.returnTimeDistribution(x, k, pairwise=pairwise, reverseComplement=reverse_complement)
