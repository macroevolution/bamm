CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(BAMM CXX)

# Set compiler flags
SET(CXX_DEBUG_FLAGS "-g")
SET(CXX_WARNING_FLAGS "-Wall -Wextra")
SET(CXX_OPTIMIZATION_FLAGS "-O3")
ADD_DEFINITIONS(${CXX_DEBUG_FLAGS}
                ${CXX_WARNING_FLAGS}
                ${CXX_OPTIMIZATION_FLAGS})

# Set Git commit id
SET(GIT_COMMIT_ID "$(shell git rev-parse HEAD)")
ADD_DEFINITIONS(-DGIT_COMMIT_ID=\"${GIT_COMMIT_ID}\")

ADD_EXECUTABLE(bamm
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
