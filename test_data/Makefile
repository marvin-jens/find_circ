.DELETE_ON_ERROR:

GENOME_HG19=/data/BIO2/pcp/systems/hg19/genome/hg19.fna
INDEX_HG19=/data/BIO2/pcp/systems/hg19/indices/bt2_hg19

cdr1as_test: cdr1as_test_out/splice_sites.bed
	@printf "\n>>> comparing to reference result. You should see 'overlap 1' here.\n\n"
	../cmp_bed.py cdr1as_test_out/splice_sites.bed cdr1as_reference_result.bed
	
clean:
	rm bt2*
	rm *.bam
	rm -rf test_out

bt2_cdr1as_locus.4.bt2: CDR1as_locus.fa
	@printf ">>> building bowtie2 index\n\n"
	bowtie2-build CDR1as_locus.fa bt2_cdr1as_locus > bt2_build.log 2>&1;

cdr1as_test.bam: bt2_cdr1as_locus.4.bt2 cdr1as_reads.fa
	@printf "\n>>> aligning example reads\n\n"
	bowtie2 -p8 --very-sensitive --reorder --mm --score-min=C,-15,0 \
		-f -U cdr1as_reads.fa -x bt2_cdr1as_locus \
		2> bt2_firstpass.log | samtools view -hbuS - | samtools sort - cdr1as_test

unmapped_cdr1as_test.bam: cdr1as_test.bam
	@printf "\n>>> fetching the unmapped reads\n\n"
	samtools view -hf 4 cdr1as_test.bam | samtools view -Sb - > unmapped_cdr1as_test.bam

anchors_cdr1as_test.fastq: unmapped_cdr1as_test.bam
	@printf "\n>>> splitting into anchors\n\n"
	../unmapped2anchors.py unmapped_cdr1as_test.bam > anchors_cdr1as_test.fastq

cdr1as_test_out/splice_sites.bed: anchors_cdr1as_test.fastq CDR1as_locus.fa
	@printf "\n>>> aligning anchors and piping through find_circ.py\n\n"
	mkdir -p cdr1as_test_out
	bowtie2 -q -U anchors_cdr1as_test.fastq -x bt2_cdr1as_locus --reorder --mm --very-sensitive --score-min=C,-15,0 2> bt2_secondpass.log | \
		../find_circ.py -G CDR1as_locus.fa -n test -p cdr1as_test_ \
		--stats cdr1as_test_out/sites.log \
		--reads cdr1as_test_out/spliced_reads.fa \
		> cdr1as_test_out/splice_sites.bed



hek_test: $(GENOME_hg19) hek_test_out/sites.bed
	@printf "\n>>> comparing HEK293 test results to reference result.\n\n"

hek_test_out/splice_sites.bed: $(GENOME_hg19) $(INDEX_HG19).4.bt2 hek_anchors.fastq.gz
	@printf "\n>>> aligning HEK293 anchors and piping through find_circ.py\n\n"
	mkdir -p hek_test_out
	bowtie2 -p8 --very-sensitive --reorder --mm --score-min=C,-15,0 \
		-q -U hek_anchors.fastq.gz  -x $(INDEX_HG19) 2> bt2_secondpass.log | \
			../find_circ.py --genome $(GENOME_HG19) \
				--name=hek_test --prefix=hek_test_ \
				--stats=hek_test_out/hek_test.log \
				--reads=hek_test_out/spliced_reads.fa \
				> hek_test_out/splice_sites.bed \


extended_test: $(GENOME_hg19) $(INDEX_HG19).4.bt2 ext_test_reads.fastq ext_test_out/sites.bed
	@printf "\n>>> comparing extended test results to reference result.\n\n"

ext_test.bam: $(INDEX_HG19).4.bt2 ext_test_reads.fastq
	@printf "\n>>> aligning example reads\n\n"
	bowtie2 -q -U ext_test_reads.fastq -x $(INDEX_HG19) -p16 --very-sensitive --mm --score-min=C,-15,0 2> bt2_firstpass.log | samtools view -hbuS - | samtools sort - ext_test

unmapped_ext_test.bam: ext_test.bam
	@printf "\n>>> fetching the unmapped reads\n\n"
	samtools view -hf 4 ext_test.bam | samtools view -Sb - > unmapped_ext_test.bam

anchors_ext_test.fastq: unmapped_ext_test.bam
	@printf "\n>>> splitting into anchors\n\n"
	../unmapped2anchors.py unmapped_ext_test.bam > anchors_ext_test.fastq

ext_test_out/sites.bed: anchors_ext_test.fastq $(GENOME_HG19)
	@printf "\n>>> aligning anchors and piping through find_circ.py\n\n"
	mkdir -p ext_test_out
	bowtie2 -q -U anchors_ext_test.fastq -x $(INDEX_HG19) --reorder --mm -M20 --score-min=C,-15,0 2> bt2_secondpass.log | \
		../find_circ.py -G $(GENOME_HG19) -n test -p ext_test_ -s ext_test_out/sites.log --halfuniq \
		> ext_test_out/sites.bed \
		2> ext_test_out/sites.reads

