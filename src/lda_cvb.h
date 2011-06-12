#ifndef TOYBOX_LDACVB_H
#define TOYBOX_LDACVB_H

#include <string>
#include <vector>
#include <map>

namespace toybox {

// inner forward declaration
class LDACVB;
class Phi;
class Updater;
class Alpha;

typedef std::vector<std::vector<int> > Int2dVec;

class LDACVB{
  public:
    LDACVB();
    ~LDACVB();
    void addDocument(const std::vector<std::string> &doc);
    void setHyperParams(double a0, double b0, double lambda);
    void learn(int num_t, int max_iterate);
    void GetResult(std::vector<std::vector<double> > *theta,
                   std::vector<std::map<std::string, double> > *beta) const;
    void ShowResult() const;
    void ShowCluster() const;
  private:
    std::map<std::string, int> dict_;
    Int2dVec words_;
    Phi *phi_;
    Updater *updater_;
    Alpha *alpha_;
    double a0_;
    double b0_;
    double lambda_;
    bool zero_;
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
  private:
    Phi();
    void normalize(int j, int i);
    void Init(int *num_i);
    int num_j_;
    int *num_i_;
    int num_t_;
    double ***values_;
};

class Updater{
  public:
    virtual ~Updater(){};
    virtual void update(const Int2dVec &words, const Alpha &alpha,
                        double lambda, int j, int i, Phi *phi) = 0;
};

class CVBUpdater : public Updater{
  public:
    CVBUpdater(const Int2dVec &words, const Phi &phi , int num_v);
    ~CVBUpdater();
    void update(const Int2dVec &words, const Alpha &alpha,
                double lambda, int j, int i, Phi *phi);
  private:
    CVBUpdater();
    void Init(const Int2dVec &words, const Phi &phi);
    void Pop(const Int2dVec &words, const Phi &phi, int j, int i);
    void Push(const Int2dVec &words, const Phi &phi, int j, int i);
    int num_j_;
    int *num_i_;
    int num_t_;
    int num_v_;
    double **num_jtd_;
    double **num_dtv_;
    double *num_dtd_;
    double **var_jtd_;
    double **var_dtv_;
    double *var_dtd_;
};

class CVB0Updater : public Updater{
  public:
    CVB0Updater(const Int2dVec &words, const Phi &phi , int num_v);
    ~CVB0Updater();
    void update(const Int2dVec &words, const Alpha &alpha,
                double lambda, int j, int i, Phi *phi);
  private:
    CVB0Updater();
    void Init(const Int2dVec &words, const Phi &phi);
    void Pop(const Int2dVec &words, const Phi &phi, int j, int i);
    void Push(const Int2dVec &words, const Phi &phi, int j, int i);
    int num_j_;
    int *num_i_;
    int num_t_;
    int num_v_;
    double **num_jtd_;
    double **num_dtv_;
    double *num_dtd_;
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

#endif // TOYBOX_LDACGS_H
