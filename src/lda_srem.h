#ifndef TOYBOX_LDASREM_H
#define TOYBOX_LDASREM_H

#include <string>
#include <vector>
#include <map>

namespace toybox {

// inner forward declaration
class LDASREM;
class Phi;
class Beta;
class Lambda;
class Alpha;

typedef std::vector<int> IntVec;
typedef std::vector<std::vector<int> > Int2dVec;

class LDASREM{
  public:
    LDASREM(int num_t, int max_iterate,
            double a0, double b0, double lambda);
    ~LDASREM();
    void addDocument(const std::vector<std::string> &doc);
    void GetResult(std::vector<double> *theta,
                   std::vector<std::map<std::string, double> > *beta) const;
    void ShowResult() const;
    void ShowCluster() const;
  private:
    LDASREM();
    void Init();
    void learn(const IntVec &words);
    std::map<std::string, int> dict_;
    Phi *phi_;
    Beta *beta_;
    Alpha *alpha_;
    Lambda *lambda_;
    int num_t_;
    int num_v_;
    int max_iterate_;
    double a0_;
    double b0_;
    double lambda0_;
};

class Phi{
  public:
    Phi(int num_t);
    ~Phi();
    int num_j() const { return num_j_; };
    int num_i() const { return num_i_; };
    int num_t() const { return num_t_; };
    double Get(int i, int t) const { return values_[i][t]; };
    void reset(const IntVec &words);
    void update(const Beta &beta, const Alpha &alpha, const IntVec &words);
  private:
    Phi();
    void normalize(int i);
    int num_j_;
    int num_i_;
    int num_t_;
    double **values_;
};

class Beta{
  public:
    Beta(int num_t);
    ~Beta();
    int num_t() const { return num_t_; };
    int num_v() const { return num_v_; };
    double Get(int t, int v) const { return values_[t][v]; };
    void Set(int t, int v, double val) { values_[t][v] = val; };
    void extend(const Lambda &lambda);
    void update(const Phi &phi, const Lambda &lambda, const IntVec &words);
  private:
    Beta();
    int num_t_;
    int num_v_;
    std::vector<std::vector<double> > values_;
};

class Lambda{
  public:
    Lambda(int num_t, double lambda0);
    ~Lambda();
    int num_t() const { return num_t_; };
    int num_v() const { return num_v_; };
    double lambda0() const { return lambda0_; };
    double Get(int t, int v) const { return values_[t][v]; };
    void Set(int t, int v, double val) { values_[t][v] = val; };
    void extend(int num_v);
    void update(const Phi &phi, const IntVec &words);
  private:
    Lambda();
    int num_t_;
    int num_v_;
    std::vector<std::vector<double> > values_;
    double lambda0_;
};

class Alpha{
  public:
    Alpha(int num_t, double a0, double b0);
    ~Alpha();
    void update(const Phi &phi);
    void updateTilde(const Phi &phi);
    int num_t() const { return num_t_; };
    double Get(int t) const { return values_[t]; };
    void Set(int t, double val) { values_[t] = val; };
  private:
    Alpha();
    int num_t_;
    double *tilde_a_;
    double tilde_b_;
    double *values_;
};

}

#endif // TOYBOX_LDASREM_H
