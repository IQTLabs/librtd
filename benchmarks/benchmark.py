import rtd_phylogeny
import skbio
import librtd
import datetime

for k in range(1, 6):

    librtd_time = datetime.timedelta()
    rtd_phylogeny_time = datetime.timedelta()

    kmers = rtd_phylogeny.getKmers(k, ["A", "T", "G", "C"])

    for seq in skbio.read("test.fasta", format="fasta"):
        seq = str(seq)

        librtd_start = datetime.datetime.now()
        librtd.return_time_distribution(seq, k)
        librtd_time += datetime.datetime.now() - librtd_start

        rtd_phylogeny_start = datetime.datetime.now()
        for kmer in kmers:
            rtd_phylogeny.getMeanSD(rtd_phylogeny.getRTD(seq, kmer))
        rtd_phylogeny_time += datetime.datetime.now() - rtd_phylogeny_start

    print(f"k={k}, librtd time: {librtd_time}, rtd_phylogeny time: {rtd_phylogeny_time}")
