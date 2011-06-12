#include "lda_cvb.h"

#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <cstring>

#include <algorithm>
#include <string>
#include <vector>
#include <map>

#include <boost/math/special_functions/digamma.hpp>

namespace toybox {

LDACVB::LDACVB()
  : phi_(NULL), updater_(NULL), alpha_(NULL),
    a0_(2.0), b0_(1.0), lambda_(1.0), zero_(false)
{
}

LDACVB::~LDACVB()
{
  if (phi_ != NULL)     { delete phi_; phi_ = NULL; }
  if (updater_ != NULL) { delete updater_; updater_ = NULL; }
  if (alpha_ != NULL)   { delete alpha_; alpha_ = NULL; }
}

void LDACVB::addDocument(const std::vector<std::string> &doc)
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

void LDACVB::setHyperParams(double a0, double b0, double lambda)
{
  a0_ = a0 > 1.0 ? a0 : 2.0;
  b0_ = b0 > 0.0 ? b0 : 1.0;
  lambda_ = lambda > 0.0 ? lambda : 1.0;
  return;
}

void LDACVB::learn(int num_t, int max_iterate)
{
  if (phi_ != NULL)     { delete phi_; }
  if (updater_ != NULL) { delete updater_; }
  if (alpha_ != NULL)   { delete alpha_; }

  int num_j = words_.size();
  int num_i[num_j];
  int num_v = dict_.size();
  for (int j = 0; j < num_j; ++j) {
    num_i[j] = words_[j].size();
  }

  phi_   = new Phi(num_j, num_i, num_t);
  if (zero_) { updater_ = new CVB0Updater(words_, *phi_, num_v); }
  else { updater_ = new CVBUpdater(words_, *phi_, num_v); }
  alpha_ = new Alpha(num_t, a0_, b0_);

  alpha_->update(*phi_);

  int iterate = 0;
  if (zero_) { fprintf(stderr,"#cvb0-iterations: %d", iterate); }
  else { fprintf(stderr,"#cvb-iterations: %d", iterate); }
  while (iterate < max_iterate) {
    ++iterate;
    for (int j = 0; j < num_j; ++j) {
      for (int i = 0; i < num_i[j]; ++i) {
        updater_->update(words_, *alpha_, lambda_, j, i, phi_);
      }
    }
    alpha_->update(*phi_);

    if (iterate % 100 == 0) { fprintf(stderr,"%d", iterate); }
    else if (iterate % 10 == 0) { fprintf(stderr,".");}
  }
  fprintf(stderr,"(%d)\n", iterate); 

  return;
}

void LDACVB::GetResult(std::vector<std::vector<double> > *theta,
                       std::vector<std::map<std::string, double> > *beta) const
{
  if (theta == NULL || beta == NULL) { return; }
  if (phi_ == NULL || alpha_ == NULL) { return; }

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
      tmp_map[inv_dict[v]] = lambda_;
    }
    beta->push_back(tmp_map);
  }
  for (int j = 0; j < num_j; ++j) {
    int num_i = phi_->num_i(j);
    for (int i = 0; i < num_i; ++i) {
      int v = words_[j][i];
      for (int t = 0; t < num_t; ++t) {
        (*beta)[t][inv_dict[v]] += phi_->Get(j, i, t);
      }
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

void LDACVB::ShowResult() const
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

void LDACVB::ShowCluster() const
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
        double v = static_cast<double>(rand()) / RAND_MAX;
        values_[j][i][t] = 1.0 + v;
      }
      this->normalize(j, i);
    }
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

CVBUpdater::CVBUpdater(const Int2dVec &words, const Phi &phi, int num_v)
  : num_j_(words.size()), num_t_(phi.num_t()), num_v_(num_v)
{
  Init(words, phi);
}

CVBUpdater::~CVBUpdater()
{
  for (int t = 0; t < num_t_; ++t) {
    delete num_dtv_[t]; num_dtv_[t] = NULL;
    delete var_dtv_[t]; var_dtv_[t] = NULL;
  }
  for (int j = 0; j < num_j_; ++j) {
    delete num_jtd_[j]; num_jtd_[j] = NULL;
    delete var_jtd_[j]; var_jtd_[j] = NULL;
  }
  delete num_dtd_; num_dtd_ = NULL;
  delete var_dtd_; var_dtd_ = NULL;
  delete num_i_; num_i_ = NULL;
}

void CVBUpdater::Init(const Int2dVec &words, const Phi &phi)
{
  num_i_ = new int[num_j_];
  num_jtd_ = new double*[num_j_];
  num_dtv_ = new double*[num_t_];
  num_dtd_ = new double[num_t_];
  var_jtd_ = new double*[num_j_];
  var_dtv_ = new double*[num_t_];
  var_dtd_ = new double[num_t_];
  std::fill_n(num_dtd_, num_t_, 0.0);
  std::fill_n(var_dtd_, num_t_, 0.0);

  for (int t = 0; t < num_t_; ++t) {
    num_dtv_[t] = new double[num_v_];
    var_dtv_[t] = new double[num_v_];
    std::fill_n(num_dtv_[t], num_v_, 0.0);
    std::fill_n(var_dtv_[t], num_v_, 0.0);
  }
  for (int j = 0; j < num_j_; ++j) {
    num_i_[j] = words[j].size();
    num_jtd_[j] = new double[num_t_];
    var_jtd_[j] = new double[num_t_];
    std::fill_n(num_jtd_[j], num_t_, 0.0);
    std::fill_n(var_jtd_[j], num_t_, 0.0);
    for (int i = 0; i < num_i_[j]; ++i) {
      int v = words[j][i];
      for (int t = 0; t < num_t_; ++t) {
        double tmp = phi.Get(j, i, t);
        num_jtd_[j][t] += tmp;
        num_dtv_[t][v] += tmp;
        num_dtd_[t] += tmp;
        var_jtd_[j][t] += tmp * (1.0 - tmp);
        var_dtv_[t][v] += tmp * (1.0 - tmp);
        var_dtd_[t] += tmp * (1.0 - tmp);
      }
    }
  }
  return;
}

void CVBUpdater::Pop(const Int2dVec &words, const Phi &phi, int j, int i)
{
  int v = words[j][i];
  for (int t = 0; t < num_t_; ++t) {
    double tmp = phi.Get(j, i, t);
    num_jtd_[j][t] -= tmp;
    num_dtv_[t][v] -= tmp;
    num_dtd_[t] -= tmp;
    var_jtd_[j][t] -= tmp * (1.0 - tmp);
    var_dtv_[t][v] -= tmp * (1.0 - tmp);
    var_dtd_[t] -= tmp * (1.0 - tmp);
  }
  return;
}

void CVBUpdater::Push(const Int2dVec &words, const Phi &phi, int j, int i)
{
  int v = words[j][i];
  for (int t = 0; t < num_t_; ++t) {
    double tmp = phi.Get(j, i, t);
    num_jtd_[j][t] += tmp;
    num_dtv_[t][v] += tmp;
    num_dtd_[t] += tmp;
    var_jtd_[j][t] += tmp * (1.0 - tmp);
    var_dtv_[t][v] += tmp * (1.0 - tmp);
    var_dtd_[t] += tmp * (1.0 - tmp);
  }
  return;
}

void CVBUpdater::update(const Int2dVec &words, const Alpha &alpha,
                        double lambda, int j, int i, Phi *phi)
{
  this->Pop(words, *phi, j, i);
  int v = words[j][i];
  double vlambda = num_v_ * lambda;
  double sum = 0.0;
  double tmp, tmp1, tmp2;
  double denom_jtd, denom_dtv, denom_dtd;
  for (int t = 0; t < num_t_; ++t) {
    tmp1 = (alpha.Get(t) + num_jtd_[j][t])
           * (lambda + num_dtv_[t][v])
           / (vlambda + num_dtd_[t]);

    denom_jtd = 2*(alpha.Get(t)+num_jtd_[j][t])*(alpha.Get(t)+num_jtd_[j][t]);
    denom_dtv = 2*(alpha.Get(t)+num_dtv_[t][v])*(alpha.Get(t)+num_dtv_[t][v]);
    denom_dtd = 2*(alpha.Get(t)+num_dtd_[t])*(alpha.Get(t)+num_dtd_[t]);
    tmp2 = - (var_jtd_[j][t] / denom_jtd)
           - (var_dtv_[t][v] / denom_dtv)
           + (var_dtd_[t] / denom_dtd);
    tmp = tmp1 * exp(tmp2);
    phi->Set(j, i, t, tmp);
    sum += tmp;
  }
  for (int t = 0; t < num_t_; ++t) {
    tmp = phi->Get(j, i, t) / sum;
    phi->Set(j, i, t, tmp);
  }
  this->Push(words, *phi, j, i);
  return;
}


CVB0Updater::CVB0Updater(const Int2dVec &words, const Phi &phi, int num_v)
  : num_j_(words.size()), num_t_(phi.num_t()), num_v_(num_v)
{
  Init(words, phi);
}

CVB0Updater::~CVB0Updater()
{
  for (int t = 0; t < num_t_; ++t) {
    delete num_dtv_[t]; num_dtv_[t] = NULL;
  }
  for (int j = 0; j < num_j_; ++j) {
    delete num_jtd_[j]; num_jtd_[j] = NULL;
  }
  delete num_dtd_; num_dtd_ = NULL;
  delete num_i_; num_i_ = NULL;
}

void CVB0Updater::Init(const Int2dVec &words, const Phi &phi)
{
  num_i_ = new int[num_j_];
  num_jtd_ = new double*[num_j_];
  num_dtv_ = new double*[num_t_];
  num_dtd_ = new double[num_t_];
  std::fill_n(num_dtd_, num_t_, 0.0);

  for (int t = 0; t < num_t_; ++t) {
    num_dtv_[t] = new double[num_v_];
    std::fill_n(num_dtv_[t], num_v_, 0.0);
  }
  for (int j = 0; j < num_j_; ++j) {
    num_i_[j] = words[j].size();
    num_jtd_[j] = new double[num_t_];
    std::fill_n(num_jtd_[j], num_t_, 0.0);
    for (int i = 0; i < num_i_[j]; ++i) {
      int v = words[j][i];
      for (int t = 0; t < num_t_; ++t) {
        double tmp = phi.Get(j, i, t);
        num_jtd_[j][t] += tmp;
        num_dtv_[t][v] += tmp;
        num_dtd_[t] += tmp;
      }
    }
  }
  return;
}

void CVB0Updater::Pop(const Int2dVec &words, const Phi &phi, int j, int i)
{
  int v = words[j][i];
  for (int t = 0; t < num_t_; ++t) {
    double tmp = phi.Get(j, i, t);
    num_jtd_[j][t] -= tmp;
    num_dtv_[t][v] -= tmp;
    num_dtd_[t] -= tmp;
  }
  return;
}

void CVB0Updater::Push(const Int2dVec &words, const Phi &phi, int j, int i)
{
  int v = words[j][i];
  for (int t = 0; t < num_t_; ++t) {
    double tmp = phi.Get(j, i, t);
    num_jtd_[j][t] += tmp;
    num_dtv_[t][v] += tmp;
    num_dtd_[t] += tmp;
  }
  return;
}

void CVB0Updater::update(const Int2dVec &words, const Alpha &alpha,
                     double lambda, int j, int i, Phi *phi)
{
  this->Pop(words, *phi, j, i);
  int v = words[j][i];
  double vlambda = num_v_ * lambda;
  double sum = 0.0;
  double tmp;
  for (int t = 0; t < num_t_; ++t) {
    tmp = (alpha.Get(t) + num_jtd_[j][t])
          * (lambda + num_dtv_[t][v])
          / (vlambda + num_dtd_[t]);
    phi->Set(j, i, t, tmp);
    sum += tmp;
  }
  for (int t = 0; t < num_t_; ++t) {
    tmp = phi->Get(j, i, t) / sum;
    phi->Set(j, i, t, tmp);
  }
  this->Push(words, *phi, j, i);
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
