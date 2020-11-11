
import terminal
import os
import times
import progress
import nimpy
import strutils
import strformat
import librtd
import tables

const version = "0.0.3.1"

# Some short templates for prettier output
template styledWrite(color: ForegroundColor, messageType: string, message: string ) = 
  stderr.styledWrite(color, messageType, ": ", fgDefault, message, "\n") 
template exit(message: string) =
  styledWrite(fgRed, "Error", message)
  quit(1)
template warn(message) = styledWrite(fgYellow, "Warning", message)
template info(message) = styledWrite(fgCyan, "Info", message) 

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

proc main*(k: int, input: string, output: string = "stdout", reverseComplement: bool = false, pairwise: bool = false) {.exportpy.} =

  # check that args are valid
  if k <= 0:
    exit(&"k value must be greater than 0, got {$k}")
  if not fileExists(input):
    exit(&"File {input} does not exist")
  if output != "stdout" and fileExists(output):
    exit(&"Output file {output} already exists. To preserve data integrity, aborting.")

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
    if reverseComplement and line.count({'U', 'u'}) > 0:
      exit("Reverse complement RTD is not currently supported for RNA sequences")

  # Add additional information
  info(&"Using librtd v{version} by Benjamin D. Lee. (c) 2020 IQT Labs, LLC.")
  when not defined(release):
    warn(&"Compiled on {CompileDate} with Nim v{NimVersion} as a debug build. This will likely be slow!")
  when not defined(danger):
    warn("Not compiled as a dangerous release build. Recompile with -d:danger for maximum performance.")

  # decide whether to write to stdout or to a file depending on the args
  var f: File
  if output != "stdout":
    f = open(output, fmWrite)
  else:
    warn("Writing data to stdout")
    f = stdout
  
  # warn the user if computing pairwise RTD
  if pairwise:
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
    let rtd = $returnTimeDistribution(sequence, k, pairwise = pairwise, reverseComplement = reverseComplement)

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