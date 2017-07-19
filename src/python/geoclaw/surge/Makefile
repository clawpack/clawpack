DAT_FILES=$(patsubst books/%.txt, %.dat, $(TXT_FILES))
PNG_FILES=$(patsubst books/%.txt, %.png, $(TXT_FILES)) 

# Generate summary table. 
results.txt : $(DAT_FILES) $(ZIPF_SRC) 
	$(ZIPF_EXE) $(DAT_FILES) > $@

# Count words
.PHONY : dats 
%.dat : books/%.txt $(COUNT_SRC) 
	$(COUNT_EXE) $< $@

dats : $(DAT_FILES) 

# Make pngs 
.PHONY : pngs

%.png : $(DAT_FILES) $(PLOT_SRC) 
	$(PLOT_EXE) $< $@ 

pngs : $(PNG_FILES) 

.PHONY : clean  
clean : 
	rm -f $(DAT_FILES)  
	rm -f $(PNG_FILES)
	rm -f results.txt 

.PHONY : variables
variables:
	@echo TXT_FILES: $(TXT_FILES)
	@echo DAT_FILES: $(DAT_FILES)



.PHONY : all
all: results.txt $(PNG_FILES)
