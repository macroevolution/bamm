CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(BAMM C CXX)

SET(BAMM_VERSION 2.2.0)
SET(BAMM_VERSION_DATE 2014-09-05)

# Set compiler flags
IF(NOT WIN32)
  SET(CXX_LANGUAGE_FLAGS "-std=c++11")
  SET(CXX_DEBUG_FLAGS "-g")
  SET(CXX_WARNING_FLAGS "-Wall -Wextra")
  SET(CXX_OPTIMIZATION_FLAGS "-O3")

  # Support 64-bit only (if using Clang)
  IF(APPLE AND (${CMAKE_CXX_COMPILER_ID} MATCHES Clang))
      SET(CXX_LANGUAGE_FLAGS "-std=c++11 -stdlib=libc++")
      SET(CXX_ARCHITECTURE_FLAGS "-arch x86_64")
      SET(CMAKE_EXE_LINKER_FLAGS "-stdlib=libc++")
  ENDIF(APPLE AND (${CMAKE_CXX_COMPILER_ID} MATCHES Clang))
ELSE(NOT WIN32)
  SET(CXX_WARNING_FLAGS "/W4")
ENDIF(NOT WIN32)

ADD_DEFINITIONS(${CXX_LANGUAGE_FLAGS}
                ${CXX_DEBUG_FLAGS}
                ${CXX_WARNING_FLAGS}
                ${CXX_OPTIMIZATION_FLAGS}
                ${CXX_ARCHITECTURE_FLAGS})

# Provide BAMM version to the compiler
ADD_DEFINITIONS(-DBAMM_VERSION=\"${BAMM_VERSION}\")
ADD_DEFINITIONS(-DBAMM_VERSION_DATE=\"${BAMM_VERSION_DATE}\")

# Set Git commit id
EXECUTE_PROCESS(COMMAND "git" "rev-parse" "HEAD"
    OUTPUT_VARIABLE GIT_COMMIT_ID
    OUTPUT_STRIP_TRAILING_WHITESPACE)
ADD_DEFINITIONS(-DGIT_COMMIT_ID=\"${GIT_COMMIT_ID}\")

SET(BAMM_SRC
    src/AcceptanceDataWriter.cpp
    src/BetaInitProposal.cpp
    src/BetaShiftProposal.cpp
    src/BetaTimeModeProposal.cpp
    src/BranchEvent.cpp
    src/BranchHistory.cpp
    src/ChainSwapDataWriter.cpp
    src/CommandLineProcessor.cpp
    src/EventCountLog.cpp
    src/EventDataWriter.cpp
    src/EventNumberProposal.cpp
    src/EventNumberForBranchProposal.cpp
    src/EventParameterProposal.cpp
    src/EventRateProposal.cpp
    src/FastSimulatePrior.cpp
    src/LambdaInitProposal.cpp
    src/LambdaShiftProposal.cpp
    src/LambdaTimeModeProposal.cpp
    src/Log.cpp
    src/MCMC.cpp
    src/MCMCDataWriter.cpp
    src/MetropolisCoupledMCMC.cpp
    src/MbRandom.cpp
    src/Model.cpp
    src/MoveEventProposal.cpp
    src/MuInitProposal.cpp
    src/MuShiftProposal.cpp
    src/NewickTreeReader.cpp
    src/Node.cpp
    src/NodeStateProposal.cpp
    src/Prior.cpp
    src/Proposal.cpp
    src/Random.cpp
    src/Settings.cpp
    src/SettingsParameter.cpp
    src/SpExBranchEvent.cpp
    src/SpExEventDataWriter.cpp
    src/SpExModel.cpp
    src/Stat.cpp
    src/StdOutDataWriter.cpp
    src/TimeModeProposal.cpp
    src/Tools.cpp
    src/TraitBranchEvent.cpp
    src/TraitEventDataWriter.cpp
    src/TraitModel.cpp
    src/Tree.cpp
    src/main.cpp)

ADD_EXECUTABLE(bamm ${BAMM_SRC})

# Find and link the pthread library for g++
IF(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
  FIND_PACKAGE(Threads REQUIRED)
  IF (Threads_FOUND)
      TARGET_LINK_LIBRARIES (bamm ${CMAKE_THREAD_LIBS_INIT})
  ENDIF()
ENDIF()

INSTALL(TARGETS bamm
        RUNTIME DESTINATION bin)
