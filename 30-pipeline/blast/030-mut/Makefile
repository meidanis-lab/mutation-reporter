SHELL := /usr/bin/bash

include ../020-xml/analysis.make

include ../../../CONFIG/$(ANALYSIS)/params.make 

%.mut: ../020-xml/%.xml	
	source ../../../venv/*/activate; python ./mutation.py $(MIN_ALIGN_LEN) $(MIN_PERC_IDENT) $(MIN_DEPTH) $(MIN_VAF) $< > $@ 


