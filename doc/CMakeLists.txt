find_package(Perl)
if(PERL_FOUND)
  if(USE_FIXED_N_RND_OPT)
    set( N_FIX_CMD "-C NFIX")
  else()
    set( N_FIX_CMD "")    
  endif()
  add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/supercell.1
                     COMMAND ${CMAKE_SOURCE_DIR}/scripts/latex2man/latex2man ${N_FIX_CMD}
                             ${CMAKE_CURRENT_SOURCE_DIR}/supercell_man.tex
                             ${CMAKE_CURRENT_BINARY_DIR}/supercell.1
                     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/supercell_man.tex
                             ${CMAKE_SOURCE_DIR}/scripts/latex2man/latex2man)
  add_custom_target(doc_man ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/supercell.1)
  set(man_1 ${CMAKE_CURRENT_BINARY_DIR}/supercell.1)
else(PERL_FOUND)
  message(WARNING "Perl is not found. Documentation (Unix man) will be not created")
  set(man_1 ${CMAKE_CURRENT_SOURCE_DIR}/supercell.1)
endif(PERL_FOUND)

find_package(LATEX)
if(PDFLATEX_COMPILER)
  if(USE_FIXED_N_RND_OPT)
    set( N_FIX_CMD_LT "\"\\def\\NFix{}\\input{${CMAKE_CURRENT_SOURCE_DIR}/supercell_man.tex}\"")
  else()
    set( N_FIX_CMD_LT "\"\\input{${CMAKE_CURRENT_SOURCE_DIR}/supercell_man.tex}\"")    
  endif()

  add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/supercell_man.pdf
                     COMMAND ${PDFLATEX_COMPILER} -interaction=batchmode 
                             ${N_FIX_CMD_LT}
                     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/supercell_man.tex)
                       
  add_custom_target(doc_pdf ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/supercell_man.pdf)
else()
  message(WARNING "LATEX is not found. PDF manual will be not created")
endif()


install(FILES ${man_1} DESTINATION share/man/man1)
