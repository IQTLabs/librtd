# librtd

![CI](https://github.com/Lab41/librtd/workflows/CI/badge.svg)

This project aims to make DNA and RNA _k_-mer return time distribution analysis simple, fast, and generalizable.

## What is _k_-mer return time distribution?

Consider the DNA sequences `AAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTT` and `ATATATATATATATATATATATATATATATATAT`.
Normal _k_-mer frequency based analysis methods would treat these sequences identically.
For _k_=1, the number of `A` 1-mers is precisely equal to the number of `T` 1-mers in both sequences.
However, _k_-mer return time methods ask the following question:
**For a given _k_-mer, how close (in base pairs) is it to the next occurrence of another _k_-mer (usually the same one).**

So, going back to the example above, for *k*=1, the return times the first sequence for `A` would be 1, 1, 1...  since each `A` *k*-mer is one base away from the next `A` *k*-mer. For the second sequence, it would be 2, 2, 2... since each `A` is takes two bases to become an `A` again. 

In `librtd`, we have generalized the concept of *k*-mer return time to include the distance of a *k*-mer not only to the next occurrence of itself but also to the next occurrence of another *k*-mer. In the first sequence, the return times from `A` to `T` are 17, 16, ..., 2, 1. This is useful for studying the relationship between the location of various pairs of *k*-mers, not just individual *k*-mers.

Once the *k*-mer return times have been calculated, `librtd` can automatically compute the mean and standard deviation of the return times for each *k*-mer, allowing the distribution to be easily summarized.

