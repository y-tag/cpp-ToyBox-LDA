#include "lda_irem.h"

#include <cstdlib>

#include <vector>

int main(int argc, char **argv){

  int num_j = 12;
  int num_i = 5;
  int data[][5] = {
    {0,1,2,3,4},
    {0,1,2,3,5},
    {0,1,2,3,8},
    {0,1,2,3,9},
    {4,5,6,7,10},
    {4,5,6,7,11},
    {4,5,6,7,0},
    {4,5,6,7,1},
    {8,9,10,11,3},
    {8,9,10,11,4},
    {8,9,10,11,6},
    {8,9,10,11,7},
  };
  std::vector<std::string> words;
  words.push_back("zero"); words.push_back("one"); words.push_back("two");
  words.push_back("three"); words.push_back("four"); words.push_back("five");
  words.push_back("six"); words.push_back("seven"); words.push_back("eight");
  words.push_back("nine"); words.push_back("ten"); words.push_back("eleven");

  toybox::LDAIREM *test = new toybox::LDAIREM();
  for (int j = 0; j < num_j; ++j) {
    std::vector<std::string> tmp_vector;
    for (int i = 0; i < num_i; ++i) {
      tmp_vector.push_back(words[data[j][i]]);
    }
    test->addDocument(tmp_vector);
  }

  test->setHyperParams(1.5, 1.0, 0.01);
  test->learn(3, 100);
  test->ShowCluster();
  test->ShowResult();

  delete test;

  return 0;
}

