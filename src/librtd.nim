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
## in
## If you are computing both the regular RTD and the reverse complement RTD, for example, you may wish to use the lower level functions instead.
## The higher level `returnTimeDistribution <#returnTimeDistribution,string, Positive>`_ is overloaded such that it can also be called with a table mapping *k*-mers to their return times.
## To compute these return times with maximum performance, first compute the *k*-mer indices using `kmerIndices <#kmerIndices,string,Positive>`_:
## 
## ```nim
## let indices = kmerIndices("ATGACGA...", k)
## ```
## Then, compute the return times using functions such as `sameKmerReturnTimes <#sameKmerReturnTimes,Table[string,seq[T][int]]>`_, `pairwiseKmerReturnTimes <#pairwiseKmerReturnTimes,Table[string,seq[T][int]]>`_ and `reverseComplementReturnTimes <#reverseComplementReturnTimes,Table[string,seq[T][int]]>`_:
## 
## ```nim
## let returnTimes = sameKmerReturnTimes(indices)
## ```
## 
## Then, to compute the metrics for the return time distribution, you can call `returnTimeDistribution <#returnTimeDistribution,Table[string,seq[T][int]]>`_ on the resultant table:
## 
## ```nim
## let rtd = returnTimeDistribution(returnTimes)
## ```
## 
## ### A Note on Degenerate Bases
##
## Degenerate bases (*i.e.* non AUTGC) bases don't have a clearly defined return time distribution since they are ambiguous.
## Therefore, to prevent invalid seqeunces resulting in invalid results, the `kmers <#kmers.i,string,Positive>`_ iterator will raise a `DegenerateBaseError <#DegenerateBaseError>`_ when confronted with a degenerate sequence. 

import strutils
import strformat
import tables
import stats

const version = "0.1"

type
  InvalidKmerLengthError* = object of CatchableError ## \
  ## Raised when the *k* value passed is too large for the given sequence.
  DegenerateBaseError* = object of CatchableError ## \
  ## Raised when an input sequence is degenerate.

iterator kmers*(x: string, k: Positive, degeneratesAllowed = false): (int, string) =
  ## Yields all of the *k*-mers (and their indices) in a given string.
  ## Note that all yielded *k*-mers are uppercase, regardless of whether the input sequence is uppercase.
  ## 
  ## This iterator can raise `InvalidKmerLengthError <#InvalidKmerLengthError>`_ when `k > x.len` as well as `DegenerateBaseError <#DegenerateBaseError>`_ when there is a non AUTGC base in the inout sequence.
  ## To override the `DegenerateBaseError <#DegenerateBaseError>`_ in the case where you explicitly want degenerate bases, call this function with `degeneratesAllowed = true`.
  runnableExamples:
    var example = newSeq[string]()
    for i, kmer in kmers("ATgCCaGA", 2):
      example.add(kmer)
    assert example == @["AT", "TG", "GC", "CC", "CA", "AG", "GA"]

  # check to make sure the *k* value isn't too big
  if k > x.len:
    raise newException(InvalidKmerLengthError, &"Unable to generate {k}-mers since {k} is longer than the input sequence, which is {x.len} bases long")

  if x.toUpper.count({'A'..'Z', '0'..'9'} - {'A', 'U', 'T', 'G',  'C'}) > 0 and not degeneratesAllowed:
    raise newException(DegenerateBaseError, "Degenerate bases do not have defined RTD.")

  for i in 0..(x.len - k):
    yield (i, x[i ..< i + k].toUpper)

func kmerIndices*(x: string, k: Positive): Table[string, seq[int]] =
  ## Returns a Table mapping *k*-mers to their indices in the input string.
  ## 
  ## In the when a *k*-mer is not present within the input string, it **will not** be in the resultant table.
  ## This is useful to prevent very large sparse tables as the value of *k* increases.
  runnableExamples:
    import tables
    let dna = "ATCGGGACCT"
    assert kmerIndices(dna, 1) == {"C": @[2, 7, 8], "A": @[0, 6], "T": @[1, 9], "G": @[3, 4, 5]}.toTable

  for i, kmer in kmers(x, k):
    if result.hasKeyOrPut(kmer, @[i]):
      result[kmer].add(i)

func sameKmerReturnTimes*(indices: Table[string, seq[int]]): Table[string, seq[int]] =
  ## Computes the return time distances (in bases) for `indices`.
  runnableExamples:
    import tables
    assert sameKmerReturnTimes(kmerIndices("ATCACA", 1)) == {"A": @[3, 2], "C": @[2]}.toTable
    
  for kmer, kmerIndices in indices.pairs:
    # there are no return times if the k-mer shows up once
    if kmerIndices.len <= 1:
      continue
    
    for i, kmerIndex in kmerIndices:
      # the last time a k-mer shows up, it doesn't have a return time
      if i == kmerIndices.high:
        break

      # otherwise, insert the return time
      if result.hasKeyOrPut(kmer, @[kmerIndices[i+1] - kmerIndices[i]]):
        result[kmer].add(kmerIndices[i+1] - kmerIndices[i])

func sameKmerReturnTimes*(x: string, k: Positive): Table[string, seq[int]] =
  ## The same function as above, but overloaded to automatically call `kmerIndices <#kmerIndices,string,Positive>`_.
  ##
  ## This is less efficient when reanalyzing the same sequence since the *k*-mer indices are recomputed.
  runnableExamples:
    import tables
    assert sameKmerReturnTimes("ATCACA", 1) == {"A": @[3, 2], "C": @[2]}.toTable
    
  sameKmerReturnTimes(kmerIndices(x, k))

func distToNextGreaterIndex*(indicies1: seq[int], indices2: seq[int]): seq[int] =
  ## Given two seqs of *k*-mer indices, calculate the distance between occurrences of the first *k*-mer and the second.
  ## 
  ## In this case, `indices1` are the indices of the first *k*-mer and `indices2` are the indices of the second *k*-mer.
  ## For each index of `indices1`, the distance to the next greater index of `indices2` will be added to the resultant seq.
  ## Note that the length of the resultant integer sequence cannot immediately be inferred, as in the example below.
  ## In the example, the greatest index of `indices1`, 9, is larger than any index of `indices2`.
  ## Therefore, it has no distance to any index of `indices2`.
  ## 
  ## You probably won't need to call this directly, but it may be useful for creating other RTD metrics.
  runnableExamples:
    assert distToNextGreaterIndex(@[1, 3, 5, 6, 9], @[2, 4, 7]) == @[1, 1, 2, 1]

  # A variable that will be used to reduce the amount of searching to find the next greater index.
  var last_indices2_idx = 0

  for idx1 in indicies1:
    # here, note that we're starting from `last_indices2_idx` in `indices2`, which is faster
    for idx2 in indices2[last_indices2_idx..indices2.high]:
      if idx2 > idx1:
        result.add(idx2 - idx1)
        break
      else:
        last_indices2_idx += 1

func pairwiseKmerReturnTimes*(indices: Table[string, seq[int]]): Table[string, seq[int]] =
  ## Calculates the return times between each pair of *k*-mers in the input table.
  ## 
  ## This can be **slow**! 
  ## Computing pairwise distances is exponential in the number of *k*-mers, which is in turn exponential in the value of *k*.
  ## 
  ## The keys in the output table are of the form `"{kmer1}_{kmer2}"` to ensure consistency with the other `KmerReturnTimes` functions, which use strings as keys.
  runnableExamples:
    import tables
    assert pairwiseKmerReturnTimes(kmerIndices("ATAAT", 1)) == {"A_T": @[1, 2, 1], "A_A": @[2, 1], "T_A": @[1], "T_T": @[3]}.toTable

  for kmer1, value1 in indices.pairs:
    for kmer2, value2 in indices.pairs:
      var distances = distToNextGreaterIndex(value1, value2)
      if distances.len > 0:
        result[&"{kmer1}_{kmer2}"] = distances

func pairwiseKmerReturnTimes*(x: string, k: Positive): Table[string, seq[int]] = 
  ## The same function as above, but overloaded to automatically call `pairwiseKmerReturnTimes <#pairwiseKmerReturnTimes,Table[string,seq[T][int]]>`_ on the result of `kmerIndices <#kmerIndices,string,Positive>`_.
  ## 
  ## As with the other overloaded functions, this may be slower when reanaylzing the same sequence since the *k*-mer indices are recomputed.
  runnableExamples:
    import tables
    assert pairwiseKmerReturnTimes("ATAAT", 1) == {"A_T": @[1, 2, 1], "A_A": @[2, 1], "T_A": @[1], "T_T": @[3]}.toTable

  pairwiseKmerReturnTimes(kmerIndices(x, k))

func reverseComplement(seq: string): string =
  ## Computes the reverse complement of a nondegenerate sequence.
  const mapping = {'A': 'T', 'a': 'T',
                   'T': 'A', 't': 'A',
                   'G': 'C', 'g': 'C',
                   'C': 'G', 'c': 'G'}.toTable
  for i in countdown(seq.high, seq.low):
    result.add(mapping[seq[i]])

func reverseComplementReturnTimes*(indices: Table[string, seq[int]]): Table[string, seq[int]] =
  ## Computes the distance from a *k*-mer to its reverse complement given a mapping of *k*-mers to their indices.
  ## 
  ## Note that that this is **not currently defined** for RNA sequences.
  runnableExamples:
    import tables
    assert reverseComplementReturnTimes(kmerIndices("ATATCCGG", 2)) == {"AT": @[2], "CC": @[2]}.toTable 

  for kmer, value in indices.pairs:
    if not indices.hasKey(kmer.reverseComplement):
      continue
    var distances = distToNextGreaterIndex(indices[kmer], indices[kmer.reverseComplement])
    if distances.len > 0:
      result[kmer] = distances

func reverseComplementReturnTimes*(x: string, k: Positive): Table[string, seq[int]] =
  ## The same as above but overloaded to automatically call reverseComplementReturnTimes <#reverseComplementReturnTimes,Table[string,seq[T][int]]>`_ after computing `kmerIndices <#kmerIndices,string,Positive>`_.
  runnableExamples:
    import tables
    assert reverseComplementReturnTimes("ATATCCGG", 2) == {"AT": @[2], "CC": @[2]}.toTable 

  reverseComplementReturnTimes(kmerIndices(x, k))

func returnTimeDistribution*(returnTimes: Table[string, seq[int]]): Table[string, float] =
  ## Given a mapping of *k*-mers to their return times, compute the mean and standard deviation of the return times.
  ## 
  ## The output table will be of the form `{"{kmer}_mean": ..., "{kmer}_std"...}` with each *k*-mer represented by two keys, one for the mean and the other for the standard deviation.
  runnableExamples:
    import tables
    assert returnTimeDistribution(sameKmerReturnTimes("AAATAGA", 1)) == {"A_mean": 1.5, "A_std": 0.5}.toTable
  for kmer, value in returnTimes.pairs:
    var statistics: RunningStat
    statistics.push(value)
    result[&"{kmer}_mean"] = statistics.mean
    result[&"{kmer}_std"] = statistics.standardDeviation

func returnTimeDistribution*(x: string, k: Positive, pairwise = false, reverseComplement = false): Table[string, float] =
  ## The master function for `librtd`, capable of accessing all of the library's functionality.
  ## 
  ## This overloaded function is capable of computing the RTD for same *k*-mers,
  ## reverse complement *k*-mers, and pairwise *k*-mers, depending on the arguments.
  ## It automatically computes summary statistics.
  ## Note that pairwise and reverse complement cannot be true.
  runnableExamples:
    import tables
    assert returnTimeDistribution("AAATAGA", 1) == {"A_mean": 1.5, "A_std": 0.5}.toTable

  if pairwise and reverseComplement:
    raise newException(ValueError, "Both pairwise and reverseComplement cannot be true")
  if pairwise:
    result = returnTimeDistribution(pairwiseKmerReturnTimes(x, k))
  elif reverseComplement:
    result = returnTimeDistribution(reverseComplementReturnTimes(x, k))
  else:
    result = returnTimeDistribution(sameKmerReturnTimes(x, k))

when isMainModule:
  import docopt
  import terminal
  import os
  import times
  import progress

  let doc = """
  Return time distribution (RTD) calculation.

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

  let args = docopt(doc, version = version)

  # parse the inputs into variables
  let input = $args["<input>"]
  let output = $args["<output>"]

  template exit(errorMesage: string) =
    ## A handy template for creating errors and exiting
    stderr.styledWrite(fgRed, "Error: ", fgDefault, errorMesage, "\n")
    quit(1)

  # attempt to parse the k value
  var k: int
  try:
    k = ($args["<k>"]).parseInt
  except ValueError:
    let invalidK = $args["<k>"]
    exit(&"k value must be an integer, not \"{invalidK}\"")

  # check that args are valid
  if k <= 0:
    exit(&"k value must be greater than 0, got {k}")
  if not fileExists(input):
    exit(&"File {input} does not exist")
  if fileExists(output):
    exit(&"Output file {output} already exists. To preserve data integrity, aborting.")

  iterator fasta(filename: string): tuple[id: string, sequence: string] =
    ## Iterate over the lines in a FASTA file, yielding one record at a time 
    var id = ""
    var row = ""
    for line in filename.lines:
      if line.startsWith(">"):
        if row != "":
          yield (id, row)
          row = ""
        id = line[1..line.high]
      else:
        row &= line.strip
    yield (id, row)
    row = ""

  # check that the sequences are non-degenerate and count how many there are
  var invalidId = ""
  var totalRecords = 0
  for line in input.lines:
    if line.startswith(">"):
      invalidId = line[1..line.high]
      totalRecords += 1
      continue
    if line.count({'a'..'z', 'A'..'Z', '0'..'9'} - {'a', 'A', 'u', 'U', 't', 'T', 'g', 'G', 'c', 'C'}) > 0:
       exit(&"Invalid (non AUTGCautgc) character in record #{totalRecords}: {invalidId}")
    
    # also check that there are no Us in the sequence if doing reverse complement RTD
    if args["--reverse-complement"] and line.count({'U', 'u'}) > 0:
      exit("Reverse complement RTD is not currently supported for RNA sequences")

  stderr.styledWrite(fgCyan, "Info: ", fgDefault, &"Using librtd v{version} by Benjamin D. Lee. (c) 2020 IQT Labs, LLC.\n")

  template warn(message) =
    stderr.styledWrite(fgYellow, "Warning: ", fgDefault, message, "\n") 
    

  # decide whether to write to stdout or to a file depending on the args
  var f: File
  if args["<output>"]:
    f = open(output, fmWrite)
  else:
    warn("Writing data to stdout")
    f = stdout
  
  # warn the user if computing pairwise RTD
  if args["--pairwise"]:
    warn("Computing pairwise RTD is much slower than computing regular or reverse complement RTD.")

  # if k > 6, warn the iser
  if k > 6:
    warn("Values of k larger than six tend to be very sparse, so use at your own peril.")

  # collect metadata about the processing
  var totalSequences = 0
  var totalBases = 0
  let time = cpuTime()

  # set up the progress bar
  var bar = newProgressBar(total=totalRecords, output=stderr)
  bar.start()

  # perform the actual computation
  for id, sequence in fasta(input):
    let rtd = $returnTimeDistribution(sequence, k, pairwise = args["--pairwise"], reverseComplement = args["--reverse-complement"])

    # if the RTD is blank, it won't be valid json, so we have to override it
    var jsonl: string
    if rtd != "{:}": 
      jsonl = "{" & &"\"id\": \"{id}\", {rtd[1..rtd.high]}" 
    else: 
      jsonl = "{" & &"\"id\": \"{id}\"" & "}" 
    f.writeLine(jsonl)

    # update the metadata
    totalSequences += 1
    totalBases += sequence.len
    bar.increment()
    
  # clean up
  f.close()
  stderr.write() # used to put the next line onto its own line from the progress bar
  stderr.styledWrite(fgGreen,
                        "Success: ",
                        fgDefault,
                        &"Analyzed RTD for {totalSequences} sequence",
                        if totalSequences > 1: "s " else: " ",
                        &"totaling {totalBases} bp ",
                        &"in {(cpuTime() - time):.1f} seconds.\n")
