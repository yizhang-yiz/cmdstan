# The bib file
file(RELATIVE_PATH BIBFILES ${CMAKE_CURRENT_SOURCE_DIR} ${Stan_SOURCE_DIR}/src/docs/bibtex/all.bib)

# All inputed texfiles
set(INPUT_TEXFILES cmdstan-guide.tex compiling.tex dump-format.tex
                   getting-started.tex installation.tex introduction.tex
                   licensing.tex output-analysis.tex references.tex running.tex
                   title.tex)

# Root tex file
set(MAIN_TEXFILE cmdstan-guide.tex)

# Files which should be configured by cmake, e.g., to add the source dir path
set(CONFIG_TEXFILES cmdstan-guide.tex title.tex references.tex)

add_latex_document(${MAIN_TEXFILE}
                   BIBFILES ${BIBFILES}
                   INPUTS ${INPUT_TEXFILES}
                   IMAGE_DIRS ${IMAGE_DIRS}
                   CONFIGURE ${CONFIG_TEXFILES}
                   USE_INDEX DEFAULT_PDF MANGLE_TARGET_NAMES)

