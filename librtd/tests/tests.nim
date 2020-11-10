# This is just an example to get you started. You may wish to put all of your
# tests into a single file, or separate them into multiple `test1`, `test2`
# etc. files (better names are recommended, just make sure the name starts with
# the letter 't').
#
# To run these tests, simply execute `nimble test`.

import unittest
import sequtils
import random
import tables
import algorithm
include librtd # needed since kmers is not exported

suite "kmers":
  test "total number of k-mers is correct":
    check: 
      toSeq(kmers("ATGCAGATA", 1)).len == 9
      toSeq(kmers("ATGCAGATA", 2)).len == 8
    var randomDNA = ""
    for i in 0..<100:
      randomDNA.add(sample(["A", "T", "G", "C"]))
    for k in 1..99:
      check:
        toSeq(kmers(randomDNA, k)).len == 100 - k + 1

  test "check the k-mers themselves are valid":
    check:
      toSeq(kmers("ATGC", 1)) == @[(0, "A"), (1, "T"), (2, "G"), (3, "C")]
      toSeq(kmers("ATGC", 2)) == @[(0, "AT"), (1, "TG"), (2, "GC")]
      toSeq(kmers("ATGC", 3)) == @[(0, "ATG"), (1, "TGC")]
      toSeq(kmers("ATGC", 4)) == @[(0, "ATGC")]

  test "automatic capitalization":
    check:
      toSeq(kmers("aTGC", 1)) == @[(0, "A"), (1, "T"), (2, "G"), (3, "C")]
      toSeq(kmers("AtGC", 2)) == @[(0, "AT"), (1, "TG"), (2, "GC")]
      toSeq(kmers("ATgC", 3)) == @[(0, "ATG"), (1, "TGC")]
      toSeq(kmers("atGc", 4)) == @[(0, "ATGC")]

  test "too large a k value raises an error":
    expect InvalidKmerLengthError:
      discard toSeq(kmers("ATGC", 10))

  test "degenerate sequence raises an error":
    expect DegenerateBaseError:
      discard toSeq(kmers("ATGCN", 1))

  test "backwards iteration":
    check:
      toSeq(kmers("ATGC", 1, fromEnd=true)) == @[(0, "A"), (1, "T"), (2, "G"), (3, "C")].reversed
      toSeq(kmers("ATGC", 2, fromEnd=true)) == @[(0, "AT"), (1, "TG"), (2, "GC")].reversed
      toSeq(kmers("ATGC", 3, fromEnd=true)) == @[(0, "ATG"), (1, "TGC")].reversed
      toSeq(kmers("ATGC", 4, fromEnd=true)) == @[(0, "ATGC")].reversed
    
suite "sameKmerReturnTimes":
  test "k=1":
    check:
      sameKmerReturnTimes("ATCACA", 1) == {"A": @[3, 2], "C": @[2]}.toTable
      sameKmerReturnTimes("CTCGTG", 1) == {"C": @[2], "T": @[3], "G": @[2]}.toTable
  test "k=2":
    check sameKmerReturnTimes("ATCGAT", 2) == {"AT": @[4]}.toTable
  test "k=3":
    check sameKmerReturnTimes("ATCATC", 3) == {"ATC": @[3]}.toTable

suite "pairwiseKmerReturnTimes":
  test "documentation example":
    check pairwiseKmerReturnTimes("ATAAT", 1) == {"A_T": @[1, 2, 1], "A_A": @[1, 2], "T_A": @[1], "T_T": @[3]}.toTable
  test "k=2 tiny example":
    check:
      pairwiseKmerReturnTimes("ATAAT", 2) == {"AT_TA": @[1], 
                                              "AT_AA": @[2], 
                                              "AT_AT": @[3],
                                              "AA_AT": @[1],
                                              "TA_AA": @[1],
                                              "TA_AT": @[2]}.toTable

suite "reverseComplementReturnTimes":
  test "documentation example":
     check reverseComplementReturnTimes("ATATCCGG", 2) == {"AT_rc": @[2], "CC_rc": @[2]}.toTable 
  test "no reverse complements in the sequence":
    for k in 2..5:
      check:
        reverseComplementReturnTimes("ATGCCCCCC", k).len == 0
  test "k=1":
    check reverseComplementReturnTimes("ATGCAGAT", 1) == {"A_rc": @[1, 3, 1],
                                                          "T_rc": @[3],
                                                          "G_rc": @[1],
                                                          "C_rc": @[2]}.toTable
  test "k=2":
    check reverseComplementReturnTimes("ATGACAGATACCACCAATGACAG", 2) == {"TG_rc": @[3, 3], 
                                                                         "AT_rc": @[9, 7], 
                                                                         "CA_rc": @[3, 6, 13]}.toTable

suite "returnTimeDistribution":
  test "documentation example":
    check returnTimeDistribution("AAATAGA", 1) == {"A_mean": 1.5, "A_std": 0.5}.toTable
  
  test "reverse complement":
    check returnTimeDistribution("ATATGGGGGGGAT", 2, reverseComplement=true) == {"AT_rc_mean": 5.5, "AT_rc_std": 3.5}.toTable

  test "check length is as expected":
    let dna = "ATAGCAAGGATACAGATA"
    for k in 2..8:
      check:
        returnTimeDistribution(dna, k).len == 2 * sameKmerReturnTimes(dna, k).len
        returnTimeDistribution(dna, k, pairwise = true).len == 2 * pairwiseKmerReturnTimes(dna, k).len
        returnTimeDistribution(dna, k, reverseComplement = true).len == 2 * reverseComplementReturnTimes(dna, k).len

    test "pairwise and reverseComplement can't both be true":
      expect ValueError:
        discard returnTimeDistribution("ATGC", 1, reverseComplement = true, pairwise = true)