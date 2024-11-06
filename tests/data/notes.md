To download the sequence for the tests

1. Download the binary from ucsc tools for your operating system, for a mac with arm processor:

```sh
wget http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.arm64/twoBitToFa
```

2. Update permissions to execute the binary:

```sh
chmod a+x twoBitToFa
```

3. Download the fasta for the region of interest:

```sh
./twoBitToFa http://hgdownload.cse.ucsc.edu/gbdb/hg19/hg19.2bit test.fa -seq=chr1 -start=250 -end=650
```
