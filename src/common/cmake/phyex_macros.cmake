macro(fortran_modules_install)
  # Récupération des arguments sous forme de liste (simples, pas nommés ici)
  set(options)
  set(one_value_args TARGET MODULE_DIR INSTALL_DIR)
  set(multi_value_args)
  cmake_parse_arguments(FMI "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  # Vérification des arguments obligatoires
  foreach(arg TARGET MODULE_DIR)
    if(NOT FMI_${arg})
      message(FATAL_ERROR "Missing required argument: ${arg}")
    endif()
  endforeach()

  # Configuration des modules Fortran
  set_target_properties(${FMI_TARGET} PROPERTIES Fortran_MODULE_DIRECTORY ${FMI_MODULE_DIR})
  target_include_directories(${FMI_TARGET} PUBLIC
    $<BUILD_INTERFACE:${FMI_MODULE_DIR}>
  )

  # Installation optionnelle
  if(FMI_INSTALL_DIR)
    target_include_directories(${FMI_TARGET} PUBLIC
      $<INSTALL_INTERFACE:${FMI_INSTALL_DIR}>
    )
    install(DIRECTORY ${FMI_MODULE_DIR}/
      DESTINATION ${FMI_INSTALL_DIR}
      COMPONENT modules
    )
  endif()
endmacro()
