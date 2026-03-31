add_executable(yvyra
    bayes.c
    charweight.c
    command.c
    diagnostics.c
    likelihood.c
    mcmc.c
    model.c
    proposal.c
    sumpt.c
    utils.c
    yamlparser.c
    yamlcmd.c
)

target_include_directories(yvyra PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}   # for generated config.h
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Math library
if(MATH_LIB)
    target_link_libraries(yvyra PRIVATE ${MATH_LIB})
endif()

# Readline
if(READLINE_LIB)
    target_link_libraries(yvyra PRIVATE ${READLINE_LIB})
endif()

# MPI
if(MPI_ENABLED)
    target_link_libraries(yvyra PRIVATE MPI::MPI_C)
endif()

install(TARGETS yvyra DESTINATION bin)
