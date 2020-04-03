dune_python_add_test(NAME pyrefelement
                     COMMAND ${PYTHON_EXECUTABLE} refelement.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_python_add_test(NAME pygeometrytype
                     COMMAND ${PYTHON_EXECUTABLE} geometrytype.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_python_add_test(NAME pytestquad
                     COMMAND ${PYTHON_EXECUTABLE} test_quad.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
