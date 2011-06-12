#ifndef TOYBOX_LDACGS_H
#define TOYBOX_LDACGS_H

#include <string>
#include <vector>
#include <map>

namespace toybox {

// inner forward declaration
class LDACGS;
class Cluster;
class Beta;
class Theta;

typedef std::vector<std::vector<int> > Int2dVec;

class LDACGS{
  public:
    LDACGS();
    ~LDACGS();
    void addDocument(const std::vector<std::string> &doc);
    void prepare(int num_t, double alpha, double lambda);
    void clear();
    void sampling();
    void GetResult(std::vector<std::vector<double> > *theta,
                   std::vector<std::map<std::string, double> > *beta) const;
    void ShowResult() const;
    void ShowCluster() const;
  private:
    std::map<std::string, int> dict_;
    Int2dVec words_;
    Cluster *cluster_;
    double alpha_;
    double lambda_;
};

class Cluster{
  public:
    Cluster(const Int2dVec &words, int num_t, int num_v);
    ~Cluster();
    int num_j() const { return num_j_; };
    int num_i(int j) const { return num_i_[j]; };
    int num_t() const { return num_t_; };
    int num_v() const { return num_v_; };
    int GetSample(const Int2dVec &words, int j, int i,
                  double alpha, double lambda) const;
    int Pop(const Int2dVec &words, int j, int i);
    void Push(const Int2dVec &words, int j, int i, int t);
    int Get(int j, int i) const { return cluster_[j][i]; };
    void Show();
  private:
    Cluster();
    void Init(const Int2dVec &words);
    int num_j_;
    int *num_i_;
    int num_t_;
    int num_v_;
    int **num_jtd_;
    int **num_dtv_;
    int *num_dtd_;
    int **cluster_;
    double *tmp_array_;
};

}

#endif // TOYBOX_LDACGS_H
