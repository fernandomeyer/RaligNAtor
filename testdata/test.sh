#! /bin/sh

set -e -x

TESTINDEX="./index/test.fa"
TEMPINDEX="./index/temp.fa"

../bin/sufconstruct ./test.fa -rna -lgslink -s $TEMPINDEX
cmp -s $TESTINDEX.alph $TEMPINDEX.alph
cmp -s $TESTINDEX.base $TEMPINDEX.base
cmp -s $TESTINDEX.des $TEMPINDEX.des
cmp -s $TESTINDEX.lcp $TEMPINDEX.lcp
cmp -s $TESTINDEX.seq $TEMPINDEX.seq
cmp -s $TESTINDEX.suf $TEMPINDEX.suf
cmp -s $TESTINDEX.sufinv $TEMPINDEX.sufinv
cmp -s $TESTINDEX.tseq $TEMPINDEX.tseq

../bin/RaligNAtorMT ./test.fa -pat ./test.pat -rna -lscan -cost 4 -indels 2 -for -rev -threads 1 -parts 1 > ./temp_online.out
diff ./result.out ./temp_online.out

../bin/RaligNAtorMT $TESTINDEX -pat ./test.pat -lgslink -cost 4 -indels 2 -for -rev -byseq -threads 1 -parts 1 > ./temp_indexed.out
diff ./result.out ./temp_indexed.out

echo "Test finished."

