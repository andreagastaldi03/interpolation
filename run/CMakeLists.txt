set(FILE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")


set(RUN_SOURCES  
   DY_example.cc
)

find_package(Ceres QUIET)
if (NOT Ceres_FOUND)
    message(WARNING "Ceres is required but was not found!")
else()
foreach(runsource ${RUN_SOURCES})
   GET_FILENAME_COMPONENT(filename ${runsource} NAME_WE)
   add_executable(${filename} ${filename}.cc ${source_files})

   target_include_directories(${filename}
      PRIVATE
      ${PROJECT_SOURCE_DIR}/src
      ${CMAKE_CURRENT_BINARY_DIR}
   )
   message(STATUS "Adding program for ${filename}")


   target_link_libraries(${filename} PUBLIC 
      m
      Interpolation::Interpolation 
      Ceres::ceres
   )

endforeach()
endif()
