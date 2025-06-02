# Simple LD R2 calculator


Docs will follow if anyone asks for them (create a github issue if so), but in the mean time:

```
git clone https://github.com/kdm9/calcLD
cd calcLD
make
./bin/calcld -s 1000000 -t $NCPUS -w 10000 -m 30 $VCF $OUTDIR
```

to make a binned version:

```
cat $OUTDIR/*.tsv | ./bin/binLD /dev/stdin $SNPWINDOW > binned.tsv
```
