#include "lda_srem.h"

#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <cstring>

#include <algorithm>
#include <string>
#include <vector>
#include <map>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

namespace toybox {

LDASREM::LDASREM(int num_t, int max_iterate,
                 double a0, double b0, double lambda)
  : phi_(NULL), beta_(NULL), alpha_(NULL), lambda_(NULL),
    num_t_(num_t), num_v_(0),
    max_iterate_(max_iterate),
    a0_(a0), b0_(b0), lambda0_(lambda)
{
  Init();
}

LDASREM::~LDASREM()
{
  if (phi_ != NULL)    { delete phi_; phi_ = NULL; }
  if (beta_ != NULL)   { delete beta_; beta_ = NULL; }
  if (alpha_ != NULL)  { delete alpha_; alpha_ = NULL; }
  if (lambda_ != NULL) { delete lambda_; lambda_ = NULL; }
}

void LDASREM::Init()
{
  phi_ = new Phi(num_t_);
  beta_ = new Beta(num_t_);
  alpha_ = new Alpha(num_t_, a0_, b0_);
  lambda_ = new Lambda(num_t_, lambda0_);
}

void LDASREM::addDocument(const std::vector<std::string> &doc)
{
  std::vector<int> words;
  
  for (std::vector<std::string>::const_iterator it = doc.begin();
       it != doc.end(); ++it) {
    int v;
    std::map<std::string, int>::iterator map_it = dict_.find(*it);
    
    if (map_it == dict_.end()) {
      v = dict_.size();
      dict_[*it] = v;
    } else {
      v = map_it->second;
    }
    words.push_back(v);
  }
  num_v_ = dict_.size();

  phi_->reset(words);
  learn(words);

  return;
}

void LDASREM::learn(const IntVec &words)
{
  int iterate = 0;
  fprintf(stderr, "%3d:", phi_->num_j());
  lambda_->extend(num_v_);
  beta_->extend(*lambda_);
  while (iterate < max_iterate_) {
    phi_->update(*beta_, *alpha_, words);
    beta_->update(*phi_, *lambda_, words);
    alpha_->update(*phi_);
    if (iterate % 100 == 0) { fprintf(stderr,"%d", iterate); }
    else if (iterate % 10 == 0) { fprintf(stderr,".");}
    ++iterate;
  }
  lambda_->update(*phi_, words);
  alpha_->updateTilde(*phi_);
  fprintf(stderr, "(%d)\n", iterate); 

  return;
}

void LDASREM::GetResult(std::vector<double> *theta,
                        std::vector<std::map<std::string, double> > *beta) const
{
  if (theta == NULL || beta == NULL) { return; }
  if (phi_ == NULL || alpha_ == NULL) { return; }

  theta->clear();
  beta->clear();

  int num_t = phi_->num_t();
  int num_v = dict_.size();
  std::map<int, std::string> inv_dict;

  for (std::map<std::string, int>::const_iterator iter = dict_.begin();
       iter != dict_.end(); ++iter) {
    inv_dict[iter->second] = iter->first;
  }

  double sum = 0.0;
  for (int t = 0; t < num_t; ++t) {
    theta->push_back(alpha_->Get(t));
    sum += alpha_->Get(t);
  }
  int num_i = phi_->num_i();
  for (int i = 0; i < num_i; ++i) {
    for (int t = 0; t < num_t; ++t) {
      (*theta)[t] += phi_->Get(i, t);
      sum += phi_->Get(i, t);
    }
  }
  for (int t = 0; t < num_t; ++t) {
    (*theta)[t] /= sum;
  }

  for (int t = 0; t < num_t; ++t) {
    std::map<std::string, double> tmp_map;
    for (int v = 0; v < num_v; ++v) {
      tmp_map[inv_dict[v]] = beta_->Get(t, v);
    }
    beta->push_back(tmp_map);
  }

  return;
}

void LDASREM::ShowResult() const
{
  std::vector<double> theta;
  std::vector<std::map<std::string, double> > beta;

  GetResult(&theta, &beta);

  int num_t = beta.size();

  fprintf(stdout, "p(z|d)\n");
  fprintf(stdout, "%3d: [", phi_->num_j());
  for (int t = 0; t < num_t; ++t) {
    fprintf(stdout, "%.10f", theta[t]);
    if (t != num_t - 1) {
      fprintf(stdout, ", ");
    }
  }
  fprintf(stdout, "]\n");

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

void LDASREM::ShowCluster() const
{
  if (phi_ == NULL) { return; }
  int num_t = phi_->num_t();
  fprintf(stdout, "%3d: [", phi_->num_j());
  int num_i = phi_->num_i();
  for (int i = 0; i < num_i; ++i) {
    double max = 0.0;
    int max_t = -1;
    for (int t = 0; t < num_t; ++t) {
      if (phi_->Get(i, t) > max) {
        max = phi_->Get(i, t);
        max_t = t;
      }
    }
    fprintf(stdout, "%d", max_t);
    if (i < num_i - 1) {
      fprintf(stdout, ", ");
    }
  }
  fprintf(stdout, "]\n");
  return;
}

Phi::Phi(int num_t)
  : num_j_(0), num_i_(0), num_t_(num_t), values_(NULL)
{
}

Phi::~Phi()
{
}

void Phi::reset(const IntVec &words)
{
  if (values_ != NULL) {
    for (int i; i < num_i_; ++i) {
      delete values_[i];
    }
    delete values_;
  }

  ++num_j_;
  num_i_ = words.size();

  values_ = new double*[num_i_];
  for (int i = 0; i < num_i_; ++i) {
    values_[i] = new double[num_t_];
    for (int t = 0; t < num_t_; ++t) {
      values_[i][t] = 1.0 + (static_cast<double>(rand()) / RAND_MAX);
    }
    this->normalize(i);
  }

  return;
}

void Phi::update(const Beta &beta, const Alpha &alpha, const IntVec &words)
{
  double phi_sum = 0.0;
  double sum = 0.0;

  for (int i = 0; i < num_i_; ++i) {
    int v = words[i];
    for (int t = 0; t < num_t_; ++t) {
      phi_sum = 0.0;
      for (int li = 0; li < num_i_; ++li) {
        phi_sum += values_[li][t];
      }

      double tmp = beta.Get(t, v)
        * exp(boost::math::digamma(alpha.Get(t) + phi_sum));
      values_[i][t] = tmp;
      sum += tmp;
    }
    for (int t = 0; t < num_t_; ++t) {
      values_[i][t] /= sum;
    }
  }


  return;
}

void Phi::normalize(int i)
{
  double sum = 0.0;
  for (int t = 0; t < num_t_; ++t) {
    sum += values_[i][t];
  }
  for (int t = 0; t < num_t_; ++t) {
    values_[i][t] /= sum;
  }
  return;
}

Beta::Beta(int num_t)
  : num_t_(num_t), num_v_(0)
{
  values_.resize(num_t);
}

Beta::~Beta()
{
}

void Beta::extend(const Lambda &lambda)
{
  int num_v = lambda.num_v();
  if (num_v == num_v_) { return; }

  for (int t = 0; t < num_t_; ++t) {
    values_[t].resize(num_v, 0.0);
    double sum = 0.0;
    for (int v = 0; v < num_v; ++v) {
      values_[t][v] = lambda.Get(t, v);
      sum += lambda.Get(t, v);
    }
    for (int v = 0; v < num_v; ++v) {
      values_[t][v] /= sum;
    }
  }
  num_v_ = num_v;

  return;
}

void Beta::update(const Phi &phi, const Lambda &lambda, const IntVec &words)
{
  for (int t = 0; t < num_t_; ++t) {
    for (int v = 0; v < num_v_; ++v) {
      values_[t][v] = lambda.Get(t, v);
    }
  }
  int num_i = phi.num_i();
  for (int i = 0; i < num_i; ++i) {
    int v = words[i];
    for (int t = 0; t < num_t_; ++t) {
      values_[t][v] += phi.Get(i, t);
    }
  }
  for (int t = 0; t < num_t_; ++t) {
    double sum = 0.0;
    for (int v = 0; v < num_v_; ++v) {
      sum += values_[t][v];
    }
    for (int v = 0; v < num_v_; ++v) {
      values_[t][v] /= sum;
    }
  }

  return;
}

Lambda::Lambda(int num_t, double lambda0)
  : num_t_(num_t), num_v_(0), lambda0_(lambda0)
{
  values_.resize(num_t);
}

Lambda::~Lambda()
{
}

void Lambda::extend(int num_v)
{
  for (int t = 0; t < num_t_; ++t) {
    values_[t].resize(num_v, lambda0_);
  }
  num_v_ = num_v;
  return;
}

void Lambda::update(const Phi &phi, const IntVec &words)
{
  int num_i = phi.num_i();
  for (int i = 0; i < num_i; ++i) {
    int v = words[i];
    for (int t = 0; t < num_t_; ++t) {
      values_[t][v] += phi.Get(i, t);
    }
  }

  return;
}


Alpha::Alpha(int num_t, double a0, double b0)
  : num_t_(num_t), tilde_b_(b0)
{
  tilde_a_ = new double[num_t];
  values_ = new double[num_t];
  for (int t = 0; t < num_t; ++t) {
    tilde_a_[t] = a0;
    double tmp = (tilde_a_[t] - 1.0) / tilde_b_;
    values_[t] = tmp > 0.0 ? tmp : 1.0E-3;
  }
}

Alpha::~Alpha()
{
  delete tilde_a_; tilde_a_ = NULL;
  delete values_; values_ = NULL;
}

void Alpha::update(const Phi &phi)
{
  double sum_alpha = 0.0;
  double n[num_t_];

  for (int t = 0; t < num_t_; ++t) {
    sum_alpha += values_[t];
    n[t] = 0.0;
  }

  int num_i = phi.num_i();
  for (int i = 0; i < num_i; ++i) {
    for (int t = 0; t < num_t_; ++t) {
      n[t] += phi.Get(i, t);
    }
  }

  double b = boost::math::digamma(sum_alpha + num_i)
             - boost::math::digamma(sum_alpha);
  tilde_b_ += b;

  for (int t = 0; t < num_t_; ++t) {
    double a_t = boost::math::digamma(values_[t] + n[t])
                 - boost::math::digamma(values_[t]);
    a_t *= values_[t];
    tilde_a_[t] += a_t;

    double tmp = (tilde_a_[t] - 1.0) / tilde_b_;
    values_[t] = tmp > 0.0 ? tmp : 1.0E-3;
  }

  return;
}

void Alpha::updateTilde(const Phi &phi)
{
  double sum_alpha = 0.0;
  double n[num_t_];

  for (int t = 0; t < num_t_; ++t) {
    sum_alpha += values_[t];
    n[t] = 0.0;
  }

  int num_i = phi.num_i();
  for (int i = 0; i < num_i; ++i) {
    for (int t = 0; t < num_t_; ++t) {
      n[t] += phi.Get(i, t);
    }
  }

  double b = boost::math::digamma(sum_alpha + num_i)
             - boost::math::digamma(sum_alpha);

  for (int t = 0; t < num_t_; ++t) {
    double a_t = boost::math::digamma(values_[t] + n[t])
                 - boost::math::digamma(values_[t]);
    a_t *= values_[t];

    double tmp = (tilde_a_[t] + a_t - 1.0) / (tilde_b_ + b);
    values_[t] = tmp > 0.0 ? tmp : 1.0E-3;
  }

  return;
}



}
