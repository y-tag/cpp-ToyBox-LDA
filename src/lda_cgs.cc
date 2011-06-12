#include "lda_cgs.h"

#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>

#include <algorithm>
#include <string>
#include <vector>
#include <map>

namespace toybox {

LDACGS::LDACGS()
  : cluster_(NULL), alpha_(1.0), lambda_(1.0)
{
}

LDACGS::~LDACGS()
{
  if (cluster_ != NULL) { delete cluster_; cluster_ = NULL; }
}

void LDACGS::addDocument(const std::vector<std::string> &doc)
{
  std::vector<int> temp_array;
  
  for (std::vector<std::string>::const_iterator it = doc.begin();
       it != doc.end(); ++it) {
    int w;
    std::map<std::string, int>::iterator map_it = dict_.find(*it);
    
    if (map_it == dict_.end()) {
      w = dict_.size();
      dict_[*it] = w;
    } else {
      w = map_it->second;
    }
    temp_array.push_back(w);
  }
  words_.push_back(temp_array);
  return;
}

void LDACGS::prepare(int num_t, double alpha, double lambda)
{
  alpha_ = alpha > 0 ? alpha : 1.0;
  lambda_ = lambda > 0 ? lambda : 1.0;
  int num_v = dict_.size();
  if (cluster_ != NULL) { delete cluster_; cluster_ = NULL; }
  cluster_ = new Cluster(words_, num_t, num_v);
  return;
}

void LDACGS::sampling()
{
  if (cluster_ == NULL) { return; }

  int t = 0;
  int num_j = cluster_->num_j();
  for (int j = 0; j < num_j; ++j) {
    int num_i = cluster_->num_i(j);
    for (int i = 0; i < num_i; ++i) {
      cluster_->Pop(words_, j, i);
      t = cluster_->GetSample(words_, j, i, alpha_, lambda_);
      cluster_->Push(words_, j, i, t);
    }
  }

  return;
}

void LDACGS::GetResult(std::vector<std::vector<double> > *theta,
                       std::vector<std::map<std::string, double> > *beta) const
{
  if (theta == NULL || beta == NULL) { return; }
  if (cluster_ == NULL) { return; }

  theta->clear();
  beta->clear();

  int num_j = cluster_->num_j();
  int num_t = cluster_->num_t();
  int num_v = dict_.size();
  std::map<int, std::string> inv_dict;

  for (std::map<std::string, int>::const_iterator iter = dict_.begin();
       iter != dict_.end(); ++iter) {
    inv_dict[iter->second] = iter->first;
  }

  (*theta).resize(num_j);
  for (int j = 0; j < num_j; ++j) {
    (*theta)[j].resize(num_t, alpha_);
    int num_i = cluster_->num_i(j);
    for (int i = 0; i < num_i; ++i) {
      int t = cluster_->Get(j, i);
      (*theta)[j][t] += 1.0;
    } 
    double sum = 0.0;
    for (int t = 0; t < num_t; ++t) {
      sum += (*theta)[j][t];
    }
    for (int t = 0; t < num_t; ++t) {
      (*theta)[j][t] /= sum;
    }
  }

  (*beta).resize(num_t);
  for (int t = 0; t < num_t; ++t) {
    for (int v = 0; v < num_v; ++v) {
      (*beta)[t][inv_dict[v]] = lambda_;
    }
  }
  for (int j = 0; j < num_j; ++j) {
    int num_i = cluster_->num_i(j);
    for (int i = 0; i < num_i; ++i) {
      int v = words_[j][i];
      int t = cluster_->Get(j, i);
      (*beta)[t][inv_dict[v]] += 1.0;
    }
  }
  for (int t = 0; t < num_t; ++t) {
    double sum = 0.0;
    for (int v = 0; v < num_v; ++v) {
      sum += (*beta)[t][inv_dict[v]];
    }
    for (int v = 0; v < num_v; ++v) {
      (*beta)[t][inv_dict[v]] /= sum;
    }
  }

  return;
}

void LDACGS::ShowResult() const
{
  std::vector<std::vector<double> > theta;
  std::vector<std::map<std::string, double> > beta;

  GetResult(&theta, &beta);

  int num_j = theta.size();
  int num_t = beta.size();

  fprintf(stdout, "p(z|d)\n");
  for (int j = 0; j < num_j; ++j) {
    fprintf(stdout, "%3d: [", j);
    for (int t = 0; t < num_t; ++t) {
      fprintf(stdout, "%.10f", theta[j][t]);
      if (t != num_t - 1) {
        fprintf(stdout, ", ");
      }
    }
    fprintf(stdout, "]\n");
  }

  fprintf(stdout, "p(v|z)\n");
  for (int t = 0; t < num_t; ++t) {
    fprintf(stdout, "%3d: [", t);
    std::map<std::string, double>::iterator iter = beta[t].begin();
    while (iter != beta[t].end()) {
      fprintf(stdout, "%s:%.10f", (iter->first).c_str(), iter->second);
      ++iter;
      if (iter != beta[t].end()) {
        fprintf(stdout, ", ");
      }
    }
    fprintf(stdout, "]\n");
  }

}

void LDACGS::ShowCluster() const
{
  if (cluster_ != NULL) {
    cluster_->Show();
  }
  return;
}


Cluster::Cluster(const Int2dVec &words, int num_t, int num_v)
  : num_j_(words.size()), num_t_(num_t), num_v_(num_v)
{
  Init(words);
}

Cluster::~Cluster()
{
  for (int t = 0; t < num_t_; ++t) {
    delete num_dtv_[t]; num_dtv_[t] = NULL;
  }
  for (int j = 0; j < num_j_; ++j) {
    delete cluster_[j]; cluster_[j] = NULL;
    delete num_jtd_[j]; num_jtd_[j] = NULL;
  }
  delete cluster_; cluster_ = NULL;
  delete num_i_; num_i_ = NULL;
  delete tmp_array_; tmp_array_ = NULL;
}

void Cluster::Init(const Int2dVec &words)
{
  num_i_ = new int[num_j_];
  cluster_ = new int*[num_j_];
  num_jtd_ = new int*[num_j_];
  num_dtv_ = new int*[num_t_];
  num_dtd_ = new int[num_t_];
  tmp_array_ = new double[num_t_];
  std::fill_n(num_dtd_, num_t_, 0);

  for (int t = 0; t < num_t_; ++t) {
    num_dtv_[t] = new int[num_v_];
    std::fill_n(num_dtv_[t], num_v_, 0);
  }
  for (int j = 0; j < num_j_; ++j) {
    num_i_[j] = words[j].size();
    cluster_[j] = new int[num_i_[j]];
    num_jtd_[j] = new int[num_t_];
    std::fill_n(num_jtd_[j], num_t_, 0);
    for (int i = 0; i < num_i_[j]; ++i) {
      int t = rand() % num_t_;
      int v = words[j][i];
      cluster_[j][i] = t;
      num_jtd_[j][t] += 1;
      num_dtv_[t][v] += 1;
      num_dtd_[t] += 1;
    }
  }
  return;
}

int Cluster::GetSample(const Int2dVec &words, int j, int i,
                       double alpha, double lambda) const
{
  int v = words[j][i];
  int t = 0;
  double vlambda = num_v_ * lambda;
  double tmp = 0.0;
  double rand_value = 0.0;
  for (t = 0; t < num_t_; ++t) {
    tmp = (alpha + num_jtd_[j][t]) * (lambda + num_dtv_[t][v]) / (vlambda + num_dtd_[t]);
    tmp_array_[t] = t == 0 ? tmp : tmp + tmp_array_[t-1];
  }
  rand_value = (static_cast<double>(rand()) / RAND_MAX) * tmp_array_[num_t_-1];
  for (t = 0; t < num_t_; ++t) {
    if (rand_value < tmp_array_[t]) {
      break;
    }
  }
  return t;
}

int Cluster::Pop(const Int2dVec &words, int j, int i)
{
  int t = cluster_[j][i];
  int v = words[j][i];
  num_jtd_[j][t] -= 1;
  num_dtv_[t][v] -= 1;
  num_dtd_[t] -= 1;
  return t;
}

void Cluster::Push(const Int2dVec &words, int j, int i, int t)
{
  int v = words[j][i];
  cluster_[j][i] = t;
  num_jtd_[j][t] += 1;
  num_dtv_[t][v] += 1;
  num_dtd_[t] += 1;
  return;
}

void Cluster::Show()
{
  for (int j = 0; j < num_j_; ++j) {
    fprintf(stdout, "%3d: [", j);
    for (int i = 0; i < num_i_[j]; ++i) {
      fprintf(stdout, "%d", cluster_[j][i]);
      if (i < num_i_[j] - 1) {
        fprintf(stdout, ", ");
      }
    }
    fprintf(stdout, "]\n");
  }
  return;
}



}
