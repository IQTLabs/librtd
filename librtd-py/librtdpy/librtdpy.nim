import librtd
import tables
import nimpy

func return_time_distribution*(x: string, k: Positive, pairwise: bool = false, reverse_complement: bool = false): Table[string, float] {.exportpy.} = 
    librtd.return_time_distribution(x, k, pairwise=pairwise, reverse_complement=reverse_complement)