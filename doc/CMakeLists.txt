find_package(Perl)
if(PERL_FOUND)
  add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/supercell.1
                     COMMAND ${CMAKE_SOURCE_DIR}/scripts/latex2man/latex2man 
                             ${CMAKE_CURRENT_SOURCE_DIR}/supercell_man.tex
                             ${CMAKE_CURRENT_BINARY_DIR}/supercell.1)
  add_custom_target(doc_man ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/supercell.1)
  set(man_1 ${CMAKE_CURRENT_BINARY_DIR}/supercell.1)
else(PERL_FOUND)
  message(WARNING "Perl is not found. Documentation will be not created")
  set(man_1 ${CMAKE_CURRENT_SOURCE_DIR}/supercell.1)
endif(PERL_FOUND)

install(FILES ${man_1} DESTINATION share/man/man1)
