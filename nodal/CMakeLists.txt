SET(APPNAME dropLandingSim_nodal)

# CMakeLists file for Forward Landing Optimization

# Identify the cpp file(s) that were to be built
FILE(GLOB SOURCE_FILES *.h *.cpp *.c)
SET(SOURCE ${SOURCE_FILES})

LINK_LIBRARIES(dropCommon)

ADD_EXECUTABLE(${APPNAME} ${SOURCE})
SET_TARGET_PROPERTIES(${APPNAME} PROPERTIES PROJECT_LABEL "Executables - ${APPNAME}")