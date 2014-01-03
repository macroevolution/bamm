CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(BAMM CXX)

# Set compiler flags
SET(CXX_DEBUG_FLAGS "-g")
SET(CXX_WARNING_FLAGS "-Wall -Wextra")
SET(CXX_OPTIMIZATION_FLAGS "-O3")

# Support Mac OS X 32-bit and 64-bit architectures
IF(APPLE)
    SET(CXX_ARCHITECTURE_FLAGS "-arch i386 -arch x86_64")
ENDIF(APPLE)

ADD_DEFINITIONS(${CXX_DEBUG_FLAGS}
                ${CXX_WARNING_FLAGS}
                ${CXX_OPTIMIZATION_FLAGS}
                ${CXX_ARCHITECTURE_FLAGS})

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
    src/Settings.cpp
    src/SettingsParameter.cpp
    src/TraitBranchEvent.cpp
    src/TraitBranchHistory.cpp
    src/TraitMCMC.cpp
    src/TraitModel.cpp
    src/Tree.cpp
    src/main.cpp)

ADD_EXECUTABLE(bamm ${BAMM_SRC})
