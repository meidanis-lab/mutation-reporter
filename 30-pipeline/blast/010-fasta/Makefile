%.fasta : ../../../00-fastq/%.fastq.gz
	gunzip < $< | sed -n '1~4s/^@/>/p;2~4p' > $@

