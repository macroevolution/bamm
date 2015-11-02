#include "NewickTreeReader.h"

#include "Tree.h"
#include "Node.h"
#include "Log.h"

#include <iostream>
#include <string>
#include <cstdlib>


const std::string NewickTreeReader::_specialCharacters = ",:;()[]";
const std::string NewickTreeReader::_whitespaceCharacters = " \n\t";


void NewickTreeReader::read(std::istream &in, Tree &tree) const
{
  peekNextCharacter(in);
  Node *rootNode = readNode(in);
  ignoreNextCharacter(in); // semicolon
  tree.setRoot(rootNode);
}


Node *NewickTreeReader::readNode(std::istream &in) const
{
  Node *node = new Node();

  if(characterStartsBranchNode(peekNextCharacter(in)))
  {
    ignoreNextCharacter(in);

    Node* leftChild = readNode(in);
    node->setLfDesc(leftChild);
    leftChild->setAnc(node);

    if(characterStartsNextChildNode(peekNextCharacter(in)))
    {
      ignoreNextCharacter(in);

      Node* rightChild = readNode(in);
      node->setRtDesc(rightChild);
      rightChild->setAnc(node);
    }
    else
    {
        exitWithError("Tree is not bifurcating.");
    }

    ignoreNextCharacter(in);
  }

  node->setName(readName(in));
  node->setBrlen(readBranchLength(in));

  if (node->getLfDesc() == NULL && node->getRtDesc() == NULL)
  {
      node->setIsTip(true);
  }

  return node;
}


bool NewickTreeReader::characterStartsBranchNode(char c) const
{
  return c == '(';
}


bool NewickTreeReader::characterEndsBranchNode(char c) const
{
  return c == ')';
}


bool NewickTreeReader::characterStartsNextChildNode(char c) const
{
  return c == ',';
}


char NewickTreeReader::peekNextCharacter(std::istream &in) const
{
  while(characterIsWhitespace(in.peek()))
    in.ignore();
  return in.peek();
}


char NewickTreeReader::readNextCharacter(std::istream &in) const
{
  while(characterIsWhitespace(in.peek()))
    in.ignore();
  return (char)in.get();
}


void NewickTreeReader::ignoreNextCharacter(std::istream &in) const
{
  in.ignore();
}


bool NewickTreeReader::characterIsWhitespace(char c) const
{
  return _whitespaceCharacters.find(c) != std::string::npos;
}


std::string NewickTreeReader::readName(std::istream &in) const
{
  return readString(in);
}


std::string NewickTreeReader::readString(std::istream &in) const
{
  std::string str;

  while(!endOfInput(in) && characterIsNormal(peekNextCharacter(in)))
    str += readNextCharacter(in);

  return str;
}


bool NewickTreeReader::characterIsNormal(char c) const
{
  return _specialCharacters.find(c) == std::string::npos;
}


bool NewickTreeReader::endOfInput(std::istream &in) const
{
  return in.eof();
}


double NewickTreeReader::readBranchLength(std::istream &in) const
{
  if(!endOfInput(in) && characterStartsBranchLength(peekNextCharacter(in)))
  {
    ignoreNextCharacter(in);
    return readNumber(in);
  }

  return 0.0;
}


double NewickTreeReader::readNumber(std::istream &in) const
{
  return atof(readString(in).c_str());
}


bool NewickTreeReader::characterStartsBranchLength(char c) const
{
  return c == ':';
}
