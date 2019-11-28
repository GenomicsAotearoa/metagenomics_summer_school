# Bin dereplication

### Objectives

* Bin dereplication using **DAS_Tool**
* Evaluating bins using **CheckM**

---

### Bin dereplication using *DAS_Tool*

As we discussed in the previous exercise, we have now generated two sets of bins from the same single assembly. With this mock data set we can see that each tool recovered 10 bins, which is the number of genomes used to create this mock community. Since we **_shouldn't_** see more than 10 bins total, these tools have obviously recovered the same genomes. However, it is not clear which tool has done a better job of recruiting contigs to each bin - we very rarely expect to see the compelte genome recovered from these kinds of data, so while it is probably the case that while an equivalent bin is present in the **MetaBAT** and **MaxBin** outputs, they will probably be of differing quality.

**DAS_Tool** is a program designed to analyse the bins in each of our binning sets and determine where these equivalent pairs (or triplets if we use three binners) exist and return the 'best' one. **DAS_Tool** does not use the actual bins, but a set of text files that link contigs to their corresponding bins in each of the bin sets. We can produce these files using **bash**.       
         
---

### Evaluating bins using *CheckM*

