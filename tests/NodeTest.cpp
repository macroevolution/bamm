#include "gtest/gtest.h"
#include "Node.h"


TEST(NodeTest, Construction)
{
    Node* node = new Node();
    Node* leftNode = new Node();
    Node* rightNode = new Node();
    Node* ancNode = new Node();

    EXPECT_EQ("", node->getName());
    node->setName("TestNode");
    EXPECT_EQ("TestNode", node->getName());

    EXPECT_EQ(NULL, node->getLfDesc());
    node->setLfDesc(leftNode);
    EXPECT_EQ(leftNode, node->getLfDesc());

    EXPECT_EQ(NULL, node->getRtDesc());
    node->setRtDesc(rightNode);
    EXPECT_EQ(rightNode, node->getRtDesc());

    EXPECT_EQ(NULL, node->getAnc());
    node->setAnc(ancNode);
    EXPECT_EQ(ancNode, node->getAnc());

    EXPECT_EQ(0, node->getIndex());
    node->setIndex(10);
    EXPECT_EQ(10, node->getIndex());

    EXPECT_EQ(0.0, node->getTime());
    node->setTime(1.0);
    EXPECT_EQ(1.0, node->getTime());

    EXPECT_EQ(0.0, node->getBrlen());
    node->setBrlen(1.0);
    EXPECT_EQ(1.0, node->getBrlen());
    
    EXPECT_EQ(0, node->getTipDescCount());
    node->setTipDescCount(10);
    EXPECT_EQ(10, node->getTipDescCount());

    EXPECT_EQ(false, node->getExtantStatus());
    node->setExtantStatus(true);
    EXPECT_EQ(true, node->getExtantStatus());

    EXPECT_EQ(false, node->getIsTip());
    node->setIsTip(true);
    EXPECT_EQ(true, node->getIsTip());

    EXPECT_EQ(false, node->getIsConstant());
    node->setIsConstant(true);
    EXPECT_EQ(true, node->getIsConstant());

    EXPECT_EQ(false, node->getIsLivingTip());
    node->setIsLivingTip(true);
    EXPECT_EQ(true, node->getIsLivingTip());

    EXPECT_EQ(0.0, node->getMapStart());
    node->setMapStart(1.0);
    EXPECT_EQ(1.0, node->getMapStart());

    EXPECT_EQ(0.0, node->getMapEnd());
    node->setMapEnd(1.0);
    EXPECT_EQ(1.0, node->getMapEnd());

    delete ancNode;
    delete rightNode;
    delete leftNode;
    delete node;
}
