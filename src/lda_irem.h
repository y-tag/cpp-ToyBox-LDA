#ifndef TOYBOX_LDAIREM_H
#define TOYBOX_LDAIREM_H

#include <string>
#include <vector>
#include <map>

namespace toybox {

// inner forward declaration
class LDAIREM;
class Phi;
class Alpha;
class Beta;

typedef std::vector<std::vector<int> > Int2dVec;

class LDAIREM{
  public:
    LDAIREM();
    ~LDAIREM();
    void addDocument(const std::vector<std::string> &doc);
    void setHyperParams(double a0, double b0, double lambda);
    void learn(int num_t, int max_iterate);
    void GetResult(std::vector<std::vector<double> > *theta,
                   std::vector<std::map<std::string, double> > *beta) const;
    void ShowResult() const;
    void ShowCluster() const;
  private:
    double calcLowerBound() const;
    std::map<std::string, int> dict_;
    Int2dVec words_;
    Phi *phi_;
    Beta *beta_;
    Alpha *alpha_;
    double a0_;
    double b0_;
    double lambda_;
};

class Phi{
  public:
    Phi(int num_j, int *num_i, int num_t);
    ~Phi();
    int num_j() const { return num_j_; };
    int num_i(int j) const { return num_i_[j]; };
    int num_t() const { return num_t_; };
    double Get(int j, int i, int t) const { return values_[j][i][t]; };
    void Set(int j, int i, int t, double val) { values_[j][i][t] = val; };
    void update(int j, int i, const Beta &beta, const Alpha &alpha,
                const Int2dVec &words);
  private:
    Phi();
    void normalize(int j, int i);
    void Init(int *num_i);
    int num_j_;
    int *num_i_;
    int num_t_;
    double ***values_;
};

class Beta{
  public:
    Beta(int num_t, int num_v);
    ~Beta();
    int num_t() const { return num_t_; };
    int num_v() const { return num_v_; };
    double Get(int t, int v) const { return values_[t][v]; };
    void Set(int t, int v, double val) { values_[t][v] = val; };
    void update(const Phi &phi, double lambda, const Int2dVec &words);
  private:
    Beta();
    int num_t_;
    int num_v_;
    double **values_;
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

#endif // TOYBOX_LDAIREM_H
