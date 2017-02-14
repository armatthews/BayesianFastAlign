#include <vector>
#include <iostream>
#include <cmath>
#include "cpyp/m.h"

using std::abs;
using std::vector;
using std::cerr;
using std::endl;
using namespace cpyp;

struct alignment_prior {
  // i is index into target sentence
  // j is index into source sentence
  // m is length of target sentence
  // n is length of source sentence (not counting NULL)
  bool use_null;
  virtual ~alignment_prior() {}
  virtual double prob(unsigned short i, unsigned short j, unsigned short m, unsigned short n) const = 0;
  virtual double null_prob(unsigned short i, unsigned short m, unsigned short n) const = 0;
  virtual double log_likelihood(const vector<vector<unsigned short>>& alignments, const vector<vector<unsigned>>& src_corpus, const unsigned tgt_corpus_token_count) const = 0;
};

// A uniform pure IBM1-style prior on alignments
struct uniform_alignment_prior : alignment_prior {
  uniform_alignment_prior(bool use_null) {
    this->use_null = use_null;
  }

  ~uniform_alignment_prior() {}

  double prob(unsigned short i, unsigned short j, unsigned short m, unsigned short n) const override {
    return 1.0;
  }

  double null_prob(unsigned short i, unsigned short m, unsigned short n) const override {
    return use_null ? 1.0 : 0.0;
  }

  double log_likelihood(const vector<vector<unsigned short>>& alignments, const vector<vector<unsigned>>& src_corpus, const unsigned tgt_corpus_token_count) const override {
    return 0.0;
  }
};

// C. Dyer, V. Chahuneau, and N. Smith (2013)'s alignment prior that favors the diagonal
// See http://www.cs.cmu.edu/~nasmith/papers/dyer+chahuneau+smith.naacl13.pdf
struct diagonal_alignment_prior : alignment_prior {
  diagonal_alignment_prior(double initial_tension, double initial_p0, bool use_null) {
    this->tension = initial_tension;
    this->p0 = initial_p0;
    this->use_null = use_null;
  }

  ~diagonal_alignment_prior() {}

  double prob(unsigned short i, unsigned short j, unsigned short m, unsigned short n) const override {
    return prob(i, j, m, n, p0, tension);
  }

  double prob(unsigned short i, unsigned short j, unsigned short m, unsigned short n, double p0, double tension) const {
    assert (i > 0 && i <= m);
    assert (j > 0 && j <= n);
    return (1.0 - p0) * exp(-tension * abs(1.0 * i / m - 1.0 * j / n));
  }

  double null_prob(unsigned short i, unsigned short m, unsigned short n) const override {
    return use_null ? p0 : 0.0;
  }

  double log_likelihood(const vector<vector<unsigned short>>& alignments, const vector<vector<unsigned>>& src_corpus, const unsigned tgt_corpus_token_count) const override {
    return log_likelihood(alignments, src_corpus, tgt_corpus_token_count, p0, tension);
  }

  double geometric_sum(const double first, const double ratio, const unsigned length) const {
    assert (length >= 1);
    return first * (1.0 - pow(ratio, length)) / (1.0 - ratio);
  }

  double compute_z(const unsigned i, const unsigned m, const unsigned n, const double p0, const double tension) const {
    bool verbose = false;
    if (verbose) cerr << "i=" << i << ", m=" << m << ", n=" << n << endl;
    /*double derp_z = 0.0; // null_prob(i, m, n);
    for (unsigned j_prime = 1; j_prime <= n; j_prime++) {
      double p = prob(i + 1, j_prime, m, n, p0, tension);
      derp_z += p;
      if (verbose) cerr << j_prime << ": " << p << endl;
    }*/
    //return derp_z;
    //if (verbose) cerr << "Total: " << derp_z << endl;

    unsigned j_up = (unsigned)(1.0 * (i + 1) * n / m);
    unsigned j_down = j_up + 1;
    if (verbose) cerr << "j_up=" << j_up << ", j_down=" << j_down << endl;
    double ratio = exp(-tension / n);
    if (verbose) cerr << "ratio=" << ratio << endl;
    double z = 0.0;
    if (j_up >= 1) {
      double j_up_p = prob(i + 1, j_up, m, n, p0, tension);
      unsigned up_len = j_up;
      double j_up_total = geometric_sum(j_up_p, ratio, up_len);
      if (verbose) cerr << "j_up_total=" << j_up_total << endl;
      z += j_up_total;
    }
    if (j_down <= n) {
      double j_down_p = prob(i + 1, j_down, m, n, p0, tension);
      unsigned down_len = n - j_down + 1;
      double j_down_total = geometric_sum(j_down_p, ratio, down_len);
      if (verbose) cerr << "j_down_total=" << j_down_total << endl;
      z += j_down_total;
    }
    if (verbose) cerr << "Final Z: " << z << endl;
    //if (verbose) cerr << "Z: " << derp_z << " / " << z << "(diff=" << derp_z - z << ")" << endl;
    //exit(1);
    return z;
  }

  double log_likelihood(const vector<vector<unsigned short>>& alignments, const vector<vector<unsigned>>& src_corpus, const unsigned tgt_corpus_token_count, double p0, double tension) const {
    // TODO: This returns NaN if the corpus contains empty sentences.
    //cerr << "[llh] p0=" << p0 << ", tension=" << tension << endl;
    assert(src_corpus.size() == alignments.size());

    const double strength = 28.44; // tgt_corpus_token_count
    const double p0_alpha = 0.08 * strength;
    const double p0_beta = 0.92 * strength;
    const double tension_shape = 70.0;
    const double tension_rate = 0.1;

    double llh = Md::log_beta_density(p0, p0_alpha, p0_beta) +
                 Md::log_gamma_density(tension, tension_shape, tension_rate);
 
    for(unsigned s = 0; s < src_corpus.size(); ++s) {
      unsigned short n = src_corpus[s].size() - 1;
      unsigned short m = alignments[s].size();
      for(unsigned short i = 0; i < alignments[s].size(); ++i) {
        unsigned short j = alignments[s][i];
        assert (j != 0 || use_null);
        assert (j >= 0 && j <= n);
        if (j == 0) {
          llh += log(null_prob(i + 1, m, n));
        }
        else {
          double p = (1.0 - p0) * prob(i + 1, j, m, n, p0, tension);
          double Z = compute_z(i, m, n, p0, tension);
          double logZ = log(Z); 
          assert (p <= Z);
          llh += log(p) - logZ;
        }
      }
    }
    return llh;
  }

  template<typename Engine>
  void resample_hyperparameters(const vector<vector<unsigned short>>& alignments, const vector<vector<unsigned>>& src_corpus, const unsigned tgt_corpus_token_count, Engine& eng, const unsigned nloop = 5, const unsigned niterations = 10) {
    // TODO: Make this faster. Initial experiments show that this is taking ~40% of the total run time.
    for (unsigned iter = 0; iter < nloop; ++iter) {
      tension = slice_sampler1d([this, &alignments, &src_corpus, &tgt_corpus_token_count](double prop_tension) { return this->log_likelihood(alignments, src_corpus, tgt_corpus_token_count, p0, prop_tension); },
                              //tension, eng, 1.0 / 30.0,
                              //30.0, 0.0, niterations, 100*niterations);
                              tension, eng, 0.1,
                              14.0, 0.0, niterations, 100*niterations);

      // TODO: This could be sped up considerably by keeping the LLH around and only update the probability of the null alignment links
      p0 = slice_sampler1d([this, &alignments, &src_corpus, &tgt_corpus_token_count](double prop_p0) { return this->log_likelihood(alignments, src_corpus, tgt_corpus_token_count, prop_p0, tension); },
                         p0, eng, 0.0, 1.0, 0.0, niterations, 100*niterations);
    }

    cerr << "Resampled diagonal alignment parameters (p0=" << p0 << ",tension=" << tension  << ") = " << log_likelihood(alignments, src_corpus, tgt_corpus_token_count) << endl;
  }

  double tension;
  double p0;
};
