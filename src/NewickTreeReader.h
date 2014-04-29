#ifndef NEWICK_TREE_READER_H
#define NEWICK_TREE_READER_H

#include <iostream>
#include <string>

#include "TreeReader.h"

class Tree;
class Node;

class NewickTreeReader : public TreeReader
{
  public:

    virtual void read(std::istream &in, Tree &tree) const;

  private:

    Node *readNode(std::istream &in) const;

    bool characterStartsBranchNode(char c) const;
    bool characterEndsBranchNode(char c) const;
    bool characterStartsNextChildNode(char c) const;
    bool characterStartsBranchLength(char c) const;
    bool characterIsNormal(char c) const;
    bool characterIsWhitespace(char c) const;

    char peekNextCharacter(std::istream &in) const;
    char readNextCharacter(std::istream &in) const;
    void ignoreNextCharacter(std::istream &in) const;
    bool endOfInput(std::istream &in) const;

    std::string readName(std::istream &in) const;
    double readBranchLength(std::istream &in) const;

    std::string readString(std::istream &in) const;
    double readNumber(std::istream &in) const;

    static const std::string _specialCharacters;
    static const std::string _whitespaceCharacters;
};

#endif
