import librtd
import tables

proc reverseComplementReturnTimeDistribution*(x: cstring, k: Positive): cstring {.exportc.} = 
  let resTable = librtd.returnTimeDistribution($x, k, reverseComplement=true)
  return ($resTable).cstring
  
proc sameKmerReturnTimeDistribution*(x: cstring, k: Positive): cstring {.exportc.} = 
  let resTable = librtd.returnTimeDistribution($x, k)
  return ($resTable).cstring
  
proc pairwiseReturnTimeDistribution*(x: cstring, k: Positive): cstring {.exportc.} = 
  let resTable = librtd.returnTimeDistribution($x, k, pairwise=true)
  return ($resTable).cstring