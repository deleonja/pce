draft.pdf :: draft.tex conclusions.tex  draft.tex intro.tex math.tex qudits.tex decoh_set.tex generators.tex  kraus.tex  pce.tex vector_space.tex proofdiagonalization.tex 2_qubits.tex
notas_JA.pdf :: notas_JA.tex
outline.pdf :: outline.tex

pce.zip :: 
	zip pce.zip 2_qubits.tex decoh_set.tex intro.tex pce.tex quantum_classical.tex a_and_generators.tex draft.tex kraus.tex proofdiagonalization.tex qudits.tex conclusions.tex generators.tex math.tex prueba.tex vector_space.tex images/twoQ_generator_diagrams-crop.pdf references.bib draft.bbl
	
pce_pra.zip :: 
	latexpand draft.tex -o draft_expanded.tex
	make draft_expanded.pdf
	rm $@ 
	zip $@ draft_expanded.tex draft_expanded.bbl images/twoQ_generator_diagrams-crop.pdf
# 	rm draft_expanded.tex
	
draft_gitdiff.pdf ::
	git checkout $(gitone)
	latexpand draft.tex -o /tmp/expanded_one.tex
	git checkout master
	git checkout $(gittwo)
	latexpand draft.tex -o /tmp/expanded_two.tex
	latexdiff /tmp/expanded_one.tex /tmp/expanded_two.tex > pce_gitdiff.tex
	git checkout master
	pdflatex -draftmode pce_gitdiff.tex
	bibtex pce_gitdiff.aux
	pdflatex -draftmode pce_gitdiff.tex
	pdflatex pce_gitdiff.tex

-include ~/makefile
# Repositorio: 
#	url = git@github.com:deleonja/QC.git
