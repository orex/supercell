find_package(LATEX COMPONENTS PDFLATEX)
find_program(PYGMENTIZE_CMD NAMES pygmentize)
if(PDFLATEX_COMPILER)
  if(PYGMENTIZE_CMD)
    configure_file(supercell_tutorial.tex supercell_tutorial.tex COPYONLY)
    configure_file(supercell_tutorial.bbl supercell_tutorial.bbl COPYONLY)
    add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/supercell_tutorial.pdf
                       COMMAND ${CMAKE_SOURCE_DIR}/scripts/latex_compile.sh supercell_tutorial
                       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/supercell_tutorial.tex
                               ${CMAKE_CURRENT_SOURCE_DIR}/supercell_tutorial.bbl)
                       
    add_custom_target(tutorial_pdf ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/supercell_tutorial.pdf)
  else()
    message(WARNING "Pygmentize is not found. PDF tutorial will be not created")
  endif()
else()
  message(WARNING "LATEX is not found. PDF tutorial will be not created")
endif()
