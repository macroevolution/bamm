CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(BAMM CXX)

# Set compiler flags
SET(CXX_DEBUG_FLAGS "-g")
SET(CXX_WARNING_FLAGS "-Wall -Wextra")
SET(CXX_OPTIMIZATION_FLAGS "-O3")

# Support Mac OS X 10.5 and greater architectures (if using Clang)
IF(APPLE AND (${CMAKE_CXX_COMPILER_ID} MATCHES Clang))
    SET(CXX_ARCHITECTURE_FLAGS "-arch i386 -arch x86_64")
ENDIF(APPLE AND (${CMAKE_CXX_COMPILER_ID} MATCHES Clang))

# Statically link C++ libraries for Windows
IF(WIN32)
    SET(CXX_LIBRARY_FLAGS "-static-libgcc -static-libstdc++ -s")
ENDIF(WIN32)

ADD_DEFINITIONS(${CXX_DEBUG_FLAGS}
                ${CXX_WARNING_FLAGS}
                ${CXX_OPTIMIZATION_FLAGS}
                ${CXX_ARCHITECTURE_FLAGS}
                ${CXX_LIBRARY_FLAGS})

# Set Git commit id
SET(GIT_COMMIT_ID "$(shell git rev-parse HEAD)")
ADD_DEFINITIONS(-DGIT_COMMIT_ID=\"${GIT_COMMIT_ID}\")

SET(BAMM_SRC
    src/Autotune.cpp
    src/BranchEvent.cpp
    src/BranchHistory.cpp
    src/FastSimulatePrior.cpp
    src/Log.cpp
    src/MCMC.cpp
    src/MbRandom.cpp
    src/Model.cpp
    src/Node.cpp
    src/Prior.cpp
    src/Settings.cpp
    src/SettingsParameter.cpp
    src/Stat.cpp
    src/TraitBranchEvent.cpp
    src/TraitBranchHistory.cpp
    src/TraitMCMC.cpp
    src/TraitModel.cpp
    src/Tree.cpp
    src/main.cpp)

ADD_EXECUTABLE(bamm ${BAMM_SRC})

# Testing by manually adding libstdc++
#FIND_LIBRARY(STD stdc++)
#MESSAGE(${STD})
#TARGET_LINK_LIBRARIES(bamm ${STD})

INSTALL(TARGETS bamm
        RUNTIME DESTINATION bin)

SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")

# Concatenate version parts
SET(MAJOR ${CPACK_PACKAGE_VERSION_MAJOR})
SET(MINOR ${CPACK_PACKAGE_VERSION_MINOR})
SET(PATCH ${CPACK_PACKAGE_VERSION_PATCH})
SET(CPACK_PACKAGE_VERSION ${MAJOR}.${MINOR}.${PATCH})

IF(WIN32)
    # To be tested on a Windows machine
    SET(CPACK_NSIS_DISPLAY_NAME "BAMM")
    SET(CPACK_NSIS_MODIFY_PATH ON)
ENDIF(WIN32)

IF(APPLE)
    SET(CPACK_GENERATOR PackageMaker)
    SET(CPACK_PACKAGE_NAME bamm)
    SET(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION})
    SET(RESOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/setup/PackageMaker)
    SET(CPACK_RESOURCE_FILE_WELCOME ${RESOURCE_DIR}/welcome.txt)
    SET(CPACK_RESOURCE_FILE_README ${RESOURCE_DIR}/readme.txt)
    SET(CPACK_RESOURCE_FILE_LICENSE ${RESOURCE_DIR}/license.txt)
ENDIF(APPLE)

# Source package
SET(CPACK_SOURCE_PACKAGE_FILE_NAME
    ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION})
SET(CPACK_SOURCE_GENERATOR TGZ)
SET(CPACK_SOURCE_IGNORE_FILES .DS_Store .git build Rcode examples problems skeleton doc)

INCLUDE(CPack)
