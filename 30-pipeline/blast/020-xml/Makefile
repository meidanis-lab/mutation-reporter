include analysis.make

include ../../../CONFIG/$(ANALYSIS)/params.make 

%.xml : ../015-join/%.fasta	
	blastx -query $< -db ../../../CONFIG/blastdb/$(ANALYSIS) -out $@ -evalue $(MAX_E_VALUE) -seg no -outfmt 5\
	  -num_threads 7


