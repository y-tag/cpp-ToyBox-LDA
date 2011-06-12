#include "lda_vb.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <string>
#include <vector>
#include <map>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

namespace toybox {

LDAVB::LDAVB()
  : phi_(NULL), gamma_(NULL), mu_(NULL), alpha_(NULL),
    a0_(2.0), b0_(1.0), lambda_(1.0)
{
}

LDAVB::~LDAVB()
{
  if (phi_ != NULL)   { delete phi_; phi_ = NULL; }
  if (gamma_ != NULL) { delete gamma_; gamma_ = NULL; }
  if (mu_ != NULL)    { delete mu_; mu_ = NULL; }
  if (alpha_ != NULL) { delete alpha_; alpha_ = NULL; }
}

void LDAVB::addDocument(const std::vector<std::string> &doc)
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

void LDAVB::setHyperParams(double a0, double b0, double lambda)
{
  a0_ = a0 > 1.0 ? a0 : 2.0;
  b0_ = b0 > 0.0 ? b0 : 1.0;
  lambda_ = lambda > 0.0 ? lambda : 1.0;
  return;
}

void LDAVB::learn(int num_t, int max_iterate)
{
  if (phi_ != NULL)   { delete phi_; }
  if (gamma_ != NULL) { delete gamma_; }
  if (mu_ != NULL)    { delete mu_; }
  if (alpha_ != NULL) { delete alpha_; }

  int num_j = words_.size();
  int num_i[num_j];
  int num_v = dict_.size();
  for (int j = 0; j < num_j; ++j) {
    num_i[j] = words_[j].size();
  }

  phi_   = new Phi(num_j, num_i, num_t);
  gamma_ = new Gamma(num_j, num_t);
  mu_    = new Mu(num_t, num_v, lambda_);
  alpha_ = new Alpha(num_t, a0_, b0_);

  alpha_->update(*phi_);
  gamma_->update(*phi_, *alpha_);
  mu_->update(*phi_, words_);

  int iterate = 0;
  fprintf(stderr,"#vb-iterations: %d", iterate); //print

  while (iterate < max_iterate) {
    ++iterate;
    phi_->update(*gamma_, *mu_, words_);
    gamma_->update(*phi_, *alpha_);
    mu_->update(*phi_, words_);
    alpha_->update(*phi_);

    if (iterate % 100 == 0) { fprintf(stderr,"%d", iterate); }
    else if (iterate % 10 == 0) { fprintf(stderr,".");}
  }
  fprintf(stderr,"(%d)\n", iterate); 

  return;
}


double LDAVB::calcLowerBound()
{
  double lower_bound = 0.0;
  int num_j = phi_->num_j();
  int num_t = phi_->num_t();
  int num_v = mu_->num_v();

  //Eq(theta)[ln p(theta|alpha)] + H[q(theta)]
  for (int j = 0; j < num_j; ++j) {
    double sum_alpha = 0.0;
    double sum_gamma = 0.0;
    for (int t = 0; t < num_t; ++t) {
      lower_bound += boost::math::lgamma(gamma_->Get(j, t))
                   - boost::math::lgamma(alpha_->Get(t));
      lower_bound += (alpha_->Get(t) - gamma_->Get(j, t))
                     * gamma_->GetDG(j, t);
      sum_alpha += alpha_->Get(t);
      sum_gamma += gamma_->Get(j, t);
    }
    lower_bound += boost::math::lgamma(sum_alpha)
                 - boost::math::lgamma(sum_gamma);
  }
  //Eq(beta)[ln p(beta|lambda)] + H[q(beta)]
  for (int t = 0; t < num_t; ++t) {
    double sum_lambda = 0.0;
    double sum_mu = 0.0;
    for (int v = 0; v < num_v; ++v) {
      lower_bound += boost::math::lgamma(mu_->Get(t, v))
                   - boost::math::lgamma(mu_->lambda());
      lower_bound += (mu_->lambda() - mu_->Get(t, v))
                     * mu_->GetDG(t, v);
      sum_lambda += mu_->lambda();
      sum_mu += mu_->Get(t, v);
    }
    lower_bound += boost::math::lgamma(sum_lambda)
                 - boost::math::lgamma(sum_mu);
  }
  //Eq(z)q(theta)[ln p(z|theta)] + H[q(z)]
  for (int j = 0; j < num_j; ++j) {
    int num_i = phi_->num_i(j);
    for (int i = 0; i < num_i; ++i) {
      for (int t = 0; t < num_t; ++t) {
        lower_bound += phi_->Get(j, i, t)
                       * gamma_->GetDG(j, t);
        if (phi_->Get(j, i, t) != 0.0) {
          lower_bound -= phi_->Get(j, i, t) * log(phi_->Get(j, i, t));
        }
      }
    }
  }
  //Eq(z)q(beta)[ln p(w|z,beta)]
  for (int j = 0; j < num_j; ++j) {
    int num_i = phi_->num_i(j);
    for (int i = 0; i < num_i; ++i) {
      int w = words_[j][i];
      for (int t = 0; t < num_t; ++t) {
        lower_bound += phi_->Get(j, i, t)
                       * mu_->GetDG(t, w);
      }
    }
  }

  return lower_bound;
}

void LDAVB::GetResult(std::vector<std::vector<double> > *theta,
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

void LDAVB::ShowResult() const
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

void LDAVB::ShowCluster() const
{
  if (phi_ == NULL) { return; }
  
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
        double v = static_cast<double>(rand()) /RAND_MAX;
        values_[j][i][t] = 1.0 + v;
      }
      this->normalize(j, i);
    }
  }

  return;
}

void Phi::update(const Gamma &gamma, const Mu &mu, const Int2dVec &words)
{
  int num_t = mu.num_t();
  int num_v = mu.num_v();
  double sum_mu[num_t];

  for (int t = 0; t < num_t; ++t) {
    double sum = 0.0;
    for (int v = 0; v < num_v; ++v) {
      sum += mu.Get(t, v);
    }
    sum_mu[t] = sum;
  }

  int num_j  = this->num_j();

  for (int j = 0; j < num_j; ++j) {
    int num_i = this->num_i(j);
    for (int i = 0; i < num_i; ++i) {
      int word = words[j][i];
      for (int t = 0; t < num_t; ++t) {
        double val = boost::math::digamma(mu.Get(t, word))
                   + boost::math::digamma(gamma.Get(j, t))
                   - boost::math::digamma(sum_mu[t]);
        this->Set(j, i, t, exp(val));
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


Gamma::Gamma(int num_j, int num_t)
  : num_j_(num_j), num_t_(num_t)
{
  values_    = new double*[num_j];
  dg_values_ = new double*[num_j];

  for (int j = 0; j < num_j_; ++j) {
    values_[j]    = new double[num_t];
    dg_values_[j] = new double[num_t];
    double dg_sum = 0.0;
    for (int t = 0; t < num_t; ++t) {
      values_[j][t] = 1.0;
      dg_sum += values_[j][t];
    }
    dg_sum = boost::math::digamma(dg_sum);
    for (int t = 0; t < num_t; ++t) {
      dg_values_[j][t] = boost::math::digamma(values_[j][t]) - dg_sum;
    }
  }
}

Gamma::~Gamma()
{
  for (int j = 0; j < num_j_; ++j) {
    delete values_[j]; values_[j] = NULL;
    delete dg_values_[j]; dg_values_[j] = NULL;
  }
  delete values_; values_ = NULL;
  delete dg_values_; dg_values_ = NULL;
}

void Gamma::update(const Phi &phi, const Alpha &alpha)
{
  int num_j = this->num_j();
  int num_t = this->num_t();

  for (int j = 0; j < num_j; ++j) {
    double dg_sum = 0.0;
    for (int t = 0; t < num_t; ++t) {
      double sum = alpha.Get(t);
      int num_i = phi.num_i(j);

      for (int i = 0; i < num_i; ++i) {
        sum += phi.Get(j, i, t);
      }
      this->Set(j, t, sum);
      dg_sum += sum;
    }
    dg_sum = boost::math::digamma(dg_sum);
    for (int t = 0; t < num_t; ++t) {
      this->SetDG(j, t, boost::math::digamma(this->Get(j, t)) - dg_sum);
    }
  }
  
  return;
}


Mu::Mu(int num_t, int num_v, double lambda)
  : num_t_(num_t), num_v_(num_v), lambda_(lambda)
{
  values_ = new double*[num_t];
  dg_values_ = new double*[num_t];

  for (int t = 0; t < num_t; ++t) {
    values_[t] = new double[num_v];
    dg_values_[t] = new double[num_v];
    double dg_sum = 0.0;
    for (int v = 0; v < num_v; ++v) {
      values_[t][v] = 1.0;
      dg_sum += values_[t][v];
    }
    dg_sum = boost::math::digamma(dg_sum);
    for (int v = 0; v < num_v; ++v) {
      dg_values_[t][v] = boost::math::digamma(values_[t][v]) - dg_sum;
    }
  }
}

Mu::~Mu()
{
  for (int t = 0; t < num_t_; ++t) {
    delete values_[t]; values_[t] = NULL;
    delete dg_values_[t]; dg_values_[t] = NULL;
  }
  delete values_; values_ = NULL;
  delete dg_values_; dg_values_ = NULL;
}

void Mu::update(const Phi &phi, const Int2dVec words)
{
  int num_t = this->num_t();
  int num_v = this->num_v();
  int num_j = phi.num_j();

  for (int t = 0; t < num_t; ++t) {
    for (int v = 0; v < num_v; ++v) {
      this->Set(t, v, lambda_);
    }
    for (int j = 0; j < num_j; ++j) {
      int num_i = phi.num_i(j);
      for (int i = 0; i < num_i; ++i) {
        int word = words[j][i];
        this->Add(t, word, phi.Get(j, i, t));
      }
    }
    double dg_sum = 0.0;
    for (int v = 0; v < num_v; ++v) {
      dg_sum += this->Get(t, v);
    }
    dg_sum = boost::math::digamma(dg_sum);
    for (int v = 0; v < num_v; ++v) {
      this->SetDG(t, v, boost::math::digamma(this->Get(t, v)) - dg_sum);
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
