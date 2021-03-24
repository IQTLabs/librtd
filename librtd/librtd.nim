## A simple library to compute return time distances and distributions for genetic sequences.
## 
## ## Basic Usage
## 
## Using `librtd` is as easy as importing it and calling the `returnTimeDistribution <#returnTimeDistribution,string,Positive>`_ function on the desired sequence:
## 
## ```nim
## import librtd
## 
## echo returnTimeDistribution("ATGACAGTA...", k=2) 
## # => {"AA_mean": 43.342, "AA_std": 12.21213, ...}
## ```
## 
## ## Advanced Usage
## 
## `librtd` is not limited to computing the distance from a *k*-mer to its next occurrence.
## Indeed, it can also compute the distance between a *k*-mer and the next occurrence of its reverse complement.
## To do so, you simply need to call `returnTimeDistribution <#returnTimeDistribution,string,Positive>`_ with the `reverseComplement=true` argument:
## 
## ```nim
## import librtd
## 
## echo returnTimeDistribution("ATGACAGTA...", k=2, reverseComplement=true) 
## # => {"AA_mean": 12.432, "AA_std": 1.058, ...}
## ```
## 
## Similarly, `librtd` can compute the pairwise RTD metric, which is just the distance from every *k*-mer to the next occurrence of every other *k*-mer.
## For example, in the sequence `ATGATGA`, the distances from "A" to "T" would be one entry, from "A" to "C", "A" to "G", and so on.
## This is slow, but computes every possible RTD.
## To compute the pairwise RTD, just call `returnTimeDistribution <#returnTimeDistribution,string,Positive>`_ with the `pairwise=true` argument: 
## 
## ```nim
## import librtd
## 
## echo returnTimeDistribution("ATGACAGTA...", k=2, pairwise=true) 
## # => {"AA_AA_mean": 32.373, "AA_AA_std": 15.892, 
## # =>  "AA_AT_mean": 95.12, "AA_AT_std": 21.824, ...}
## ```
## 
## ## A Note on Degenerate Bases
##
## Degenerate bases (*i.e.* non AUTGC) bases don't have a clearly defined return time distribution since they are ambiguous.
## Therefore, to prevent invalid seqeunces resulting in invalid results, a `DegenerateBaseError <#DegenerateBaseError>`_ will be raised when confronted with a degenerate sequence. 
## 
## ## A Note on Naming
## 
## Nim and Python, although they share nearly identical syntax at the high level, vary signifcantly with respect to naming.
## In Python, `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ states that "function names should be lowercase, with words separated by underscores as necessary to improve readability."
## Nim *prefers* camelCase, but, importantly, doesn't distinguish between the two due to its definition of `identifier equality <https://nim-lang.github.io/Nim/manual.html#lexical-analysis-identifier-equality>`_.
## Therefore, calling `kmer_indices` in Nim is precisely the same as calling `kmerIndices`.
## Within the Nim implementation, the names are *defined* in Nim-style camel case.
## In the Python API, these names are exported only in the PEP8 style.
## Therefore, the Nim function `kmerIndices` is exported to Python only as `kmer_indices` for Python style compatibility.

import strutils
import strformat
import tables
import stats

type
  InvalidKmerLengthError* = object of CatchableError ## \
  ## Raised when the *k* value passed is too large for the given sequence.
  DegenerateBaseError* = object of CatchableError ## \
  ## Raised when an input sequence is degenerate.

iterator kmers(x: string, k: Positive, degeneratesAllowed = false, fromEnd = false): (int, string) =
  ## Yields all of the *k*-mers (and their indices) in a given string.
  ## 
  ## Note that all yielded *k*-mers are uppercase, regardless of whether the input sequence is uppercase.
  ## For backwards iteration, set `fromEnd=true`.
  ## 
  ## This iterator can raise `InvalidKmerLengthError <#InvalidKmerLengthError>`_ when `k > x.len` as well as `DegenerateBaseError <#DegenerateBaseError>`_ when there is a non AUTGC base in the inout sequence.
  ## To override the `DegenerateBaseError <#DegenerateBaseError>`_ in the case where you explicitly want degenerate bases, call this function with `degeneratesAllowed = true`.
  runnableExamples:
    import sequtils
    assert toSeq(kmers("ATGC", 2)) == @[(0, "AT"), (1, "TG"), (2, "GC")]
    assert toSeq(kmers("ATGC", 2, fromEnd=true)) == @[(2, "GC"), (1, "TG"), (0, "AT")]
    assert toSeq(kmers("gCCaGA", 2)).mapIt(it[1]) == @["GC", "CC", "CA", "AG", "GA"]

  # check to make sure the *k* value isn't too big
  if k > x.len:
    raise newException(InvalidKmerLengthError, &"Unable to generate {k}-mers since {k} is longer than the input sequence, which is {x.len} bases long")

  if x.toUpperAscii.count({'A'..'Z', '0'..'9'} - {'A', 'U', 'T', 'G', 'C'}) > 0 and not degeneratesAllowed:
    raise newException(DegenerateBaseError, "Degenerate bases do not have defined RTD.")

  if not fromEnd:
    for i in 0..(x.len - k):
      yield (i, x[i ..< i + k].toUpperAscii)
  else:
    for i in countdown(x.len - k, 0):
      yield (i, x[i ..< i + k].toUpperAscii) 

func sameKmerReturnTimes*(x: string, k: Positive): Table[string, seq[int]] = 
  ## Compute the return times for *k*-mers in `x`.
  ## 
  ## This function is exported to Python.
  runnableExamples:
    import tables
    assert sameKmerReturnTimes("ATCACA", 1) == {"A": @[3, 2], "C": @[2]}.toTable
  var lastIndex = initTable[string, int]()
  for i, kmer in kmers(x, k):
    if kmer in lastIndex and result.hasKeyOrPut(kmer, @[i - lastIndex[kmer]]):
        result[kmer].add(i - lastIndex[kmer])
    lastIndex[kmer] = i

proc reverseComplement*(x: cstring): cstring = 
  # Note that this function takes cstrings for compatibility with the JS backend.
  # Not sure why it needed it but it refused to work with vanilla Nim strings.
  var res = newString(x.len)
  var i = x.high
  var j = 0
  while i >= 0:
    if x[j] == 'A':
      res[i] = 'T'
    elif  x[j] == 'T':
      res[i] = 'A'
    elif  x[j] == 'G':
      res[i] = 'C'
    elif  x[j] == 'C':
      res[i] = 'G'
    else:
      raise newException(CatchableError, "Invalid character in sequence")
    inc(j)
    dec(i)
  return res.cstring

func pairwiseKmerReturnTimes*(x: string, k: Positive): Table[string, seq[int]] =
  ## Calculates the return times between each pair of *k*-mers in the input table.
  ## 
  ## This can be **slow**! 
  ## Computing pairwise distances is exponential in the number of *k*-mers, which is in turn exponential in the value of *k*.
  ## 
  ## The keys in the output table are of the form `"{kmer1}_{kmer2}"`.
  ## Note that the value sequences are in reverse order.
  ## That is, if the distances between a pair of *k*-mers are 4, 3, 5 then the value will be `@[5, 3, 4]`.
  ## This is due to backwards iteration during computation of the times and has no effect on the distribution.
  ## 
  ## This function is exported to Python. 
  runnableExamples:
    import tables
    assert pairwiseKmerReturnTimes("ATAAT", 1) == 
      {"A_T": @[1, 2, 1], "A_A": @[1, 2], "T_A": @[1], "T_T": @[3]}.toTable
  var lastIndex = initTable[string, int]()
  for i, kmer in x.kmers(k, fromEnd=true):
    for kmer2 in lastIndex.keys:
      if result.hasKeyOrPut(kmer & "_" & kmer2, @[lastIndex[kmer2] - i]):
        result[kmer & "_" & kmer2].add(lastIndex[kmer2] - i)
    lastIndex[kmer] = i

func reverseComplementReturnTimes*(x: string, k: Positive): Table[string, seq[int]] = 
  ## Computes the distance from a *k*-mer to its reverse complement given a mapping of *k*-mers to their indices.
  ## 
  ## This metric is **not currently defined for RNA sequences**.
  ## The keys in the output table are of the form `"{kmer1}_rc"`.
  ## Note that the value sequences are in reverse order.
  ## That is, if the distances between a *k*-mer and its reverse complement are 4, 3, 5 then the value will be `@[5, 3, 4]`.
  ## This is due to backwards iteration during computation of the times and has no effect on the distribution.
  ##  
  ## This function is exported to Python. 
  runnableExamples:
    import tables
    assert reverseComplementReturnTimes("ATATCCGG", 2) == 
      {"AT_rc": @[2], "CC_rc": @[2]}.toTable 
  var lastIndex = initTable[string, int]()
  var rc: string
  for i, kmer in x.kmers(k, fromEnd=true):
    rc = $(kmer.reverseComplement) # conversion needed for JS
    if rc in lastIndex and result.hasKeyOrPut(kmer & "_rc", @[lastIndex[rc] - i]):
        result[kmer & "_rc"].add(lastIndex[rc] - i)
    lastIndex[kmer] = i

func returnTimeDistribution*(returnTimes: Table[string, seq[int]]): Table[string, float] =
  ## Given a mapping of *k*-mers to their return times, compute the mean and standard deviation of the return times.
  ## 
  ## The output table will be of the form `{"{kmer}_mean": ..., "{kmer}_std"...}` with each *k*-mer represented by two keys, one for the mean and the other for the standard deviation.
  runnableExamples:
    import tables
    assert returnTimeDistribution(sameKmerReturnTimes("AAATAGA", 1))  == 
      {"A_mean": 1.5, "A_std": 0.5}.toTable
  for kmer, value in returnTimes.pairs:
    var statistics: RunningStat
    statistics.push(value)
    result[&"{kmer}_mean"] = statistics.mean
    result[&"{kmer}_std"] = statistics.standardDeviation

func returnTimeDistribution*(x: string, k: Positive, pairwise: bool = false, reverseComplement: bool = false): Table[string, float] =
  ## The master function for `librtd`, capable of accessing all of the library's functionality.
  ## 
  ## This overloaded function is capable of computing the RTD for same *k*-mers,
  ## reverse complement *k*-mers, and pairwise *k*-mers, depending on the arguments.
  ## It automatically computes summary statistics.
  ## Note that pairwise and reverse complement cannot be true.
  ## 
  ## This function is exported to Python. 
  runnableExamples:
    import tables
    assert returnTimeDistribution("AAATAGA", 1) == {"A_mean": 1.5, "A_std": 0.5}.toTable
    assert returnTimeDistribution("AAAAGACCGCC", 1, reverseComplement=true) == 
      {"G_rc_std": 0.5, "G_rc_mean": 1.5, "C_rc_mean": 1.5, "C_rc_std": 0.5}.toTable

  if pairwise and reverseComplement:
    raise newException(ValueError, "Both pairwise and reverseComplement cannot be true")
  if pairwise:
    result = returnTimeDistribution(pairwiseKmerReturnTimes(x, k))
  elif reverseComplement:
    result = returnTimeDistribution(reverseComplementReturnTimes(x, k))
  else:
    result = returnTimeDistribution(sameKmerReturnTimes(x, k))
