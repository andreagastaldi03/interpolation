set(FILE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")


set(RUN_SOURCES  
   dose_engine.cc
   dose_engine_smallStep.cc
   dose_engine_noSc.cc
)

foreach(runsource ${RUN_SOURCES})
   GET_FILENAME_COMPONENT(filename ${runsource} NAME_WE)

   # 1. dichiara eseguibile
   add_executable(${filename} ${filename}.cc)

   # 2. dichiara dove sono i file .hh
   target_include_directories(${filename}
      PRIVATE
      ${PROJECT_SOURCE_DIR}/src
      ${PROJECT_SOURCE_DIR}/include
      ${CMAKE_CURRENT_BINARY_DIR}
      /home/andreagastaldi/computational/runge-kutta/rk_cpp
   )
   message(STATUS "Adding program for ${filename}")

   # 3. collega le librerie 
   target_link_libraries(${filename} PUBLIC 
      m
      Interpolation::Interpolation 
   )

endforeach()
