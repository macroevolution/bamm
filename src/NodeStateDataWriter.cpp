#include "NodeStateDataWriter.h"
#include "Settings.h"
#include "TraitModel.h"
#include "Tree.h"

#include <iostream>


NodeStateDataWriter::NodeStateDataWriter(Settings& settings) :
    _outputFileName(settings.get("nodeStateOutfile")),
    _outputFreq(settings.get<int>("nodeStateWriteFreq"))
{
    if (_outputFreq > 0) {
        initializeStream();
    }
}


void NodeStateDataWriter::initializeStream()
{
    _outputStream.open(_outputFileName.c_str());
}


NodeStateDataWriter::~NodeStateDataWriter()
{
    if (_outputFreq > 0) {
        _outputStream.close();
    }
}


void NodeStateDataWriter::writeData(int generation, TraitModel& model)
{
    if (_outputFreq == 0 || generation % _outputFreq != 0) {
        return;
    }

    Tree& tree = *model.getTreePtr();

    _outputStream << generation;
    tree.writeBranchPhenotypes(tree.getRoot(), _outputStream);
    _outputStream << ";\n";
}
