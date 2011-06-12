#ifndef TOYBOX_LDAVB_H
#define TOYBOX_LDAVB_H

#include <string>
#include <vector>
#include <map>

namespace toybox {

// inner forward declaration
class LDAVB;
class Phi;
class Gamma;
class Mu;
class Alpha;

typedef std::vector<std::vector<int> > Int2dVec;


class LDAVB{
  public:
    LDAVB();
    ~LDAVB();
    void addDocument(const std::vector<std::string> &doc);
    void setHyperParams(double a0, double b0, double lambda);
    void learn(int num_t, int max_iterate);
    void GetResult(std::vector<std::vector<double> > *theta,
                   std::vector<std::map<std::string, double> > *beta) const;
    void ShowResult() const;
    void ShowCluster() const;
  private:
    double calcLowerBound();
    Int2dVec words_;
    std::map<std::string, int> dict_;
    Phi *phi_;
    Gamma *gamma_;
    Mu *mu_;
    Alpha *alpha_;
    double a0_;
    double b0_;
    double lambda_;
};

// q(z|phi) q(theta|gamma) q(beta|mu)

class Phi{
  public:
    Phi(int num_j, int *num_i, int num_t);
    ~Phi();
    void update(const Gamma &gamma, const Mu &mu, const Int2dVec &words);
    int num_j() const { return num_j_; };
    int num_i(int j) const { return num_i_[j]; };
    int num_t() const { return num_t_; };
    double Get(int j, int i, int t) const { return values_[j][i][t]; };
    void Set(int j, int i, int t, double val) { values_[j][i][t] = val; };
  private:
    Phi();
    void Init(int *num_i);
    void normalize(int j, int i);
    int num_j_;
    int *num_i_;
    int num_t_;
    double ***values_;
};

class Gamma{
  public:
    Gamma(int num_j, int num_t);
    ~Gamma();
    void update(const Phi &phi, const Alpha &alpha);
    int num_j() const { return num_j_; };
    int num_t() const { return num_t_; };
    double Get(int j, int t) const { return values_[j][t]; };
    void Set(int j, int t, double val) { values_[j][t] = val; };
    double GetDG(int j, int t) const { return dg_values_[j][t]; };
    void SetDG(int j, int t, double val) { dg_values_[j][t] = val; };
  private:
    Gamma();
    int num_j_;
    int num_t_;
    double **values_;
    double **dg_values_;
};

class Mu{
  public:
    Mu(int num_t, int num_v, double lambda);
    ~Mu();
    void update(const Phi &phi, const Int2dVec words);
    int num_t() const { return num_t_; };
    int num_v() const { return num_v_; };
    double lambda() const { return lambda_; };
    double Get(int t, int v) const { return values_[t][v]; };
    void Set(int t, int v, double val) { values_[t][v] = val; };
    void Add(int t, int v, double val) { values_[t][v] += val; };
    double GetDG(int t, int v) const { return dg_values_[t][v]; };
    void SetDG(int t, int v, double val) { dg_values_[t][v] = val; };
  private:
    Mu();
    int num_t_;
    int num_v_;
    double lambda_;
    double **values_;
    double **dg_values_;
};

class Alpha{
  public:
    Alpha(int num_t, double a0, double b0);
    ~Alpha();
    void update(const Phi &phi);
    int num_t() const { return num_t_; };
    double a0() const { return a0_; };
    double b0() const { return b0_; };
    double Get(int t) const { return values_[t]; };
    void Set(int t, double val) { values_[t] = val; };
  private:
    Alpha();
    int num_t_;
    double a0_;
    double b0_;
    double *values_;
};

}

#endif //TOYBOX_LDAVB_H

