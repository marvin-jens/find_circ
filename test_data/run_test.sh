#!/usr/bin/env bash
rm bt2*
rm *.bam
rm -rf test_out

echo ">>> building bowtie2 index..."
bowtie2-build CDR1as_locus.fa bt2_cdr1as_locus &> bt2_build.log
echo ">>> aligning example reads"
bowtie2 -p16 --very-sensitive --mm -M20 --score-min=C,-15,0 -x bt2_cdr1as_locus -f -U reads.fa 2> bt2_firstpass.log | samtools view -hbuS - | samtools sort - test_vs_cdr1as
echo ">>> get the unmapped"
samtools view -hf 4 test_vs_cdr1as.bam | samtools view -Sb - > unmapped_test.bam
echo ">>> split into anchors"
../unmapped2anchors.py unmapped_test.bam > test_anchors.qfa
echo ">>> run find_circ.py"
mkdir test_out
bowtie2 --reorder --mm -M20 --score-min=C,-15,0 -q -x bt2_cdr1as_locus -U test_anchors.qfa 2> bt2_secondpass.log | ../find_circ.py -r ../samples_example.txt -G . -p cdr1as_test_ -s test_out/sites.log > test_out/sites.bed 2> test_out/sites.reads
echo ">>> compare to reference result. You should see 'overlap 1' here."
../cmp_bed.py test_out/sites.bed result_test.bed
