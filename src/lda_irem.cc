#include "lda_irem.h"

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

LDAIREM::LDAIREM()
  : phi_(NULL), beta_(NULL), alpha_(NULL),
    a0_(2.0), b0_(1.0), lambda_(1.0)
{
}

LDAIREM::~LDAIREM()
{
  if (phi_ != NULL)   { delete phi_; phi_ = NULL; }
  if (beta_ != NULL)  { delete beta_; beta_ = NULL; }
  if (alpha_ != NULL) { delete alpha_; alpha_ = NULL; }
}

void LDAIREM::addDocument(const std::vector<std::string> &doc)
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

void LDAIREM::setHyperParams(double a0, double b0, double lambda)
{
  a0_ = a0 > 1.0 ? a0 : 2.0;
  b0_ = b0 > 0.0 ? b0 : 1.0;
  lambda_ = lambda > 0.0 ? lambda : 1.0;
  return;
}

void LDAIREM::learn(int num_t, int max_iterate)
{
  if (phi_ != NULL)   { delete phi_; }
  if (beta_ != NULL)  { delete beta_; }
  if (alpha_ != NULL) { delete alpha_; }

  int num_j = words_.size();
  int num_i[num_j];
  int num_v = dict_.size();
  for (int j = 0; j < num_j; ++j) {
    num_i[j] = words_[j].size();
  }

  phi_   = new Phi(num_j, num_i, num_t);
  beta_  = new Beta(num_t, num_v);
  alpha_ = new Alpha(num_t, a0_, b0_);

  beta_->update(*phi_, lambda_, words_);
  alpha_->update(*phi_);

  int iterate = 0;
  fprintf(stderr,"#irem-iterations: %d", iterate); //print

  while (iterate < max_iterate) {
    ++iterate;
    for (int j = 0; j < num_j; ++j) {
      for (int i = 0; i < num_i[j]; ++i) {
        phi_->update(j, i, *beta_, *alpha_, words_);
      }
      beta_->update(*phi_, lambda_, words_);
    }
    alpha_->update(*phi_);

    if (iterate % 100 == 0) { fprintf(stderr,"%d", iterate); }
    else if (iterate % 10 == 0) { fprintf(stderr,".");}
  }
  fprintf(stderr,"(%d)\n", iterate); 

  return;
}

double LDAIREM::calcLowerBound() const
{
  double lower_bound = 0.0;
  int num_j = phi_->num_j();
  int num_t = phi_->num_t();

  for (int j = 0; j < num_j; ++j) {
    for (int i = 0; i < phi_->num_i(j); ++j) {
      int v = words_[j][i];
      for (int t = 0; t < num_t; ++t) {
        if (phi_->Get(j, i, t) != 0.0) {
          lower_bound += phi_->Get(j, i, t)
                         * (log(beta_->Get(t, v)) - log(phi_->Get(j, i, t)));
        }
      }
    }
  }

  double sum_alpha = 0.0;
  for (int t = 0; t < num_t; ++t) {
    sum_alpha += alpha_->Get(t);
  }

  for (int j = 0; j < num_j; ++j) {
    int num_i = phi_->num_i(j);
    lower_bound += boost::math::lgamma(sum_alpha);
    lower_bound -= boost::math::lgamma(num_i + sum_alpha);
    for (int t = 0; t < num_t; ++t) {
      lower_bound -= boost::math::lgamma(alpha_->Get(t));
      double sum_phi = 0.0;
      for (int i = 0; i < phi_->num_i(j); ++i) {
        sum_phi += phi_->Get(j, i, t);
      }
      lower_bound += boost::math::lgamma(alpha_->Get(t) + sum_phi);
    }
  }

  return lower_bound;
}

void LDAIREM::GetResult(std::vector<std::vector<double> > *theta,
                        std::vector<std::map<std::string, double> > *beta) const
{
  if (theta == NULL || beta == NULL) { return; }
  if (phi_ == NULL || beta_ == NULL || alpha_ == NULL) { return; }

  theta->clear();
  beta->clear();

  int num_j = phi_->num_j();
  int num_t = phi_->num_t();
  int num_v = dict_.size();
  std::map<int, std::string> inv_dict;

  for (std::map<std::string, int>::const_iterator iter = dict_.begin();
       iter != dict_.end(); ++iter) {
    inv_dict[iter->second] = iter->first;
  }

  for (int j = 0; j < num_j; ++j) {
    std::vector<double> tmp_vec;
    double sum = 0.0;
    for (int t = 0; t < num_t; ++t) {
      tmp_vec.push_back(alpha_->Get(t));
      sum += alpha_->Get(t);
    }
    int num_i = phi_->num_i(j);
    for (int i = 0; i < num_i; ++i) {
      for (int t = 0; t < num_t; ++t) {
        tmp_vec[t] += phi_->Get(j, i, t);
        sum += phi_->Get(j, i, t);
      }
    }
    for (int t = 0; t < num_t; ++t) {
      tmp_vec[t] /= sum;
    }
    theta->push_back(tmp_vec);
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

void LDAIREM::ShowResult() const
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

void LDAIREM::ShowCluster() const
{
  if (phi_ != NULL) {
    int num_j = phi_->num_j();
    int num_t = phi_->num_t();
    for (int j = 0; j < num_j; ++j) {
      fprintf(stdout, "%3d: [", j);
      int num_i = phi_->num_i(j);
      for (int i = 0; i < num_i; ++i) {
        double max = 0.0;
        int max_t = -1;
        for (int t = 0; t < num_t; ++t) {
          if (phi_->Get(j, i, t) > max) {
            max = phi_->Get(j, i, t);
            max_t = t;
          }
        }
        fprintf(stdout, "%d", max_t);
        if (i < num_i - 1) {
          fprintf(stdout, ", ");
        }
      }
      fprintf(stdout, "]\n");
    }
  }
  return;
}

Phi::Phi(int num_j, int *num_i, int num_t)
  : num_j_(num_j), num_t_(num_t)
{
  Init(num_i);
}

Phi::~Phi()
{
  for (int j = 0; j < num_j_; ++j) {
    for (int i = 0; i < num_i_[j]; ++i) {
      delete values_[j][i]; values_[j][i] = NULL;
    }
    delete values_[j]; values_[j] = NULL;
  }
  delete values_; values_ = NULL;
}

void Phi::Init(int *num_i)
{
  num_i_ = new int[num_j_];
  values_ = new double**[num_j_];

  memmove(num_i_, num_i, sizeof(int) * num_j_);
  
  for (int j = 0; j < num_j_; ++j) {
    values_[j] = new double*[num_i[j]];
    for (int i = 0; i < num_i_[j]; ++i) {
      values_[j][i] = new double[num_t_];
      for (int t = 0; t < num_t_; ++t) {
        values_[j][i][t] = 1.0 + static_cast<double>(rand()) / RAND_MAX;
      }
      this->normalize(j, i);
    }
  }

  return;
}

void Phi::update(int j, int i, const Beta &beta, const Alpha &alpha,
                 const Int2dVec &words)
{
  int v = words[j][i];
  double phi_sum = 0.0;
  double sum = 0.0;

  for (int t = 0; t < num_t_; ++t) {
    phi_sum = 0.0;
    for (int li = 0; li < num_i_[j]; ++li) {
      phi_sum += this->Get(j, li, t);
    }

    double tmp = beta.Get(t, v)
                 * exp(boost::math::digamma(alpha.Get(t) + phi_sum));
    this->Set(j, i, t, tmp);
    sum += tmp;
  }
  for (int t = 0; t < num_t_; ++t) {
    double tmp = this->Get(j, i, t) / sum;
    this->Set(j, i, t, tmp);
  }

  return;
}

void Phi::normalize(int j, int i)
{
  double sum = 0.0;
  
  sum = 0.0;
  for (int t = 0; t < num_t_; ++t) {
    sum += values_[j][i][t];
  }
  for (int t = 0; t < num_t_; ++t) {
    values_[j][i][t] /= sum;
  }
  return;
}

Beta::Beta(int num_t, int num_v)
  : num_t_(num_t), num_v_(num_v)
{
  values_ = new double*[num_t];
  for (int t = 0; t < num_t; ++t) {
    values_[t] = new double[num_v];
    for (int v = 0; v < num_v; ++v) {
      values_[t][v] = 1.0 + static_cast<double>(rand()) / RAND_MAX;
    }
  }
}

Beta::~Beta()
{
  for (int t = 0; t < num_t_; ++t) {
    delete values_[t]; values_[t] = NULL;
  }
  delete values_; values_ = NULL;
}

void Beta::update(const Phi &phi, double lambda, const Int2dVec &words)
{
  for (int t = 0; t < num_t_; ++t) {
    std::fill_n(values_[t], num_v_, lambda);
  }
  for (int j = 0; j < phi.num_j(); ++j) {
    for (int i = 0; i < phi.num_i(j); ++i) {
      int v = words[j][i];
      for (int t = 0; t < num_t_; ++t) {
        values_[t][v] += phi.Get(j, i, t);
      }
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


Alpha::Alpha(int num_t, double a0, double b0)
  : num_t_(num_t), a0_(a0), b0_(b0)
{
  values_ = new double[num_t];
  for (int t = 0; t < num_t; ++t) {
    values_[t] = 1.0;
  }
}

Alpha::~Alpha()
{
  delete values_; values_ = NULL;
}

void Alpha::update(const Phi &phi)
{
  int num_t = this->num_t();
  int num_j = phi.num_j();

  double sum_alpha = 0.0;
  double n[num_j][num_t];

  for (int t = 0; t < num_t; ++t) {
    sum_alpha += this->Get(t);
    for (int j = 0; j < num_j; ++j) {
      n[j][t] = 0.0;
    }
  }

  for (int j = 0; j < num_j; ++j) {
    int num_i = phi.num_i(j);
    for (int i = 0; i < num_i; ++i) {
      for (int t = 0; t < num_t; ++t) {
        n[j][t] += phi.Get(j, i, t);
      }
    }
  }

  for (int t = 0; t < num_t; ++t) {
    double nume = 0.0;
    double deno = 0.0;

    double alpha_t = this->Get(t);

    for (int j = 0; j < num_j; ++j) {
      nume += boost::math::digamma(alpha_t + n[j][t])
            - boost::math::digamma(alpha_t);

      deno += boost::math::digamma(phi.num_i(j) + sum_alpha)
            - boost::math::digamma(sum_alpha);
    }

    nume = a0_ - 1 + nume * alpha_t;
    deno = b0_ + deno;

    this->Set(t, nume / deno);
  }

  return;
}

}
