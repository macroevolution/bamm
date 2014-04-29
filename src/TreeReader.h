#ifndef TREE_READER_H
#define TREE_READER_H

#include <iostream>

class Tree;

class TreeReader
{
  public:

    virtual void read(std::istream &in, Tree &tree) const = 0;
};

#endif
