#include <execinfo.h>
#include <signal.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cstdlib>
#include "boost/archive/text_oarchive.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/program_options.hpp"

#include "bilingual_corpus.h"
#include "cpyp/m.h"
#include "cpyp/random.h"
#include "cpyp/crp.h"
#include "cpyp/tied_parameter_resampler.h"
#include "alignment_prior.h"

using namespace std;
using namespace cpyp;
namespace po=boost::program_options;

double log_likelihood(const tied_parameter_resampler<crp<unsigned>>& base_ttable_params,
                      const tied_parameter_resampler<crp<unsigned>>& underlying_ttable_params, 
                      const diagonal_alignment_prior& ap,
                      const vector<vector<unsigned>>& src_corpus,
                      const unsigned tgt_corpus_token_count,
                      const crp<unsigned>& base_ttable,
                      const vector<crp<unsigned>>& underlying_ttable,
                      const vector<vector<unsigned short>>& alignments) {
  double llh = 0.0;
  llh += base_ttable_params.log_likelihood();
  llh += underlying_ttable_params.log_likelihood();

  llh += base_ttable.log_likelihood();
  for (auto& crp : underlying_ttable)
    llh += crp.log_likelihood();

  llh += ap.log_likelihood(alignments, src_corpus, tgt_corpus_token_count);
  return llh;
}

void output_alignments(vector<vector<unsigned>>& tgt_corpus, vector<vector<unsigned short>>& alignments, bool reverse) {
  for (unsigned i = 0; i < tgt_corpus.size(); i++) {
    for (unsigned j = 0; j < tgt_corpus[i].size(); j++) {
      if (alignments[i][j] != 0) {
        if (!reverse) {
          cout << alignments[i][j] - 1 << "-" << j << " ";
        }
        else {
          cout << j << "-" << alignments[i][j] - 1 << " ";
        }
      }
    }
    cout << "\n";
  }
}

void output_latent_variables(vector<crp<unsigned>>& underlying_ttable, Dict& src_dict, Dict& tgt_dict, Dict& doc_dict) {
  vector<unsigned> ind(tgt_dict.max() + 1);
  for (unsigned tgt_id = 0; tgt_id < tgt_dict.max() + 1; tgt_id++)
    ind[tgt_id] = tgt_id;

  cerr << "=====BEGIN TTABLE=====\n";
  for (unsigned src_id = 0; src_id < underlying_ttable.size(); src_id++) {
    crp<unsigned>& p = underlying_ttable[src_id];
    if (p.num_customers() == 0)
      continue;
    cerr << src_dict.Convert(src_id) << "\t" << p.num_customers()
      << "\t" << p.num_tables() << endl;
    vector<unsigned> translations;
    for (auto it = p.begin(); it != p.end(); ++it) {
      translations.push_back(it->first);
    }
    sort(translations.begin(), translations.end(), [&p](unsigned a, unsigned b) { return p.num_customers(a) > p.num_customers(b); });
    for (unsigned i = 0; i < translations.size(); ++i) {
      unsigned tgt_id = translations[i]; 
      cerr << "\t" << tgt_dict.Convert(tgt_id) << "\t" << p.prob(tgt_id, 1.0 / tgt_dict.max())
        << "\t" << p.num_customers(tgt_id) << "\t" << p.num_tables(tgt_id) << endl;
    } 
    cerr << "\t" << "[other]" << "\t" << p.prob(0, 1.0 / tgt_dict.max()) * (tgt_dict.max() - translations.size()) << endl;
  }
}

vector<vector<unsigned short>> load_alignment_file(const string& filename, const vector<vector<unsigned>>& src_corpus, const vector<vector<unsigned>>& tgt_corpus) {
  ifstream in(filename);
  string line;
  int lc = 0;
  vector<vector<unsigned short>> alignments(src_corpus.size());
  while (getline(in, line)) {
    ++lc;
    int src_size = src_corpus[lc - 1].size();
    int tgt_size = tgt_corpus[lc - 1].size();
    vector<unsigned short> alignment(tgt_size, 0);
    vector<string> links;
    boost::split(links, line, boost::is_any_of(" "));
    for (unsigned k = 0; k < links.size(); ++k) {
      int dash_location = links[k].find('-');
      int i = atoi(links[k].substr(0, dash_location).c_str());
      int j = atoi(links[k].substr(dash_location + 1).c_str());
      assert(i >= 0);
      assert(i < src_size);
      assert(j >= 0);
      assert( j < tgt_size);
      alignment[j] = i + 1;
    }
    alignments[lc - 1] = alignment;
  }
  return alignments;
}

void handler(int sig) {
        void* array[10];
        int size = backtrace(array, 10);
        cerr << "Error: signal " << sig << ":\n";
        backtrace_symbols_fd(array, size, STDERR_FILENO);
        exit(1);
}

int main(int argc, char** argv) {
  signal(SIGSEGV, handler);
  po::options_description options("Options");
  options.add_options()
    ("training_corpus,i", po::value<string>()->required(), "Training corpus, in format of source ||| target or docid ||| source ||| target")
    ("samples,n", po::value<int>()->required(), "Number of samples")
    ("alignments,a", po::value<string>(), "Initial alignments")
    ("reverse,r", "Swap source and target sides of corpus")
    ("quiet,q", "Don't output any of the latent variables -- just the alignments, please.")
    ("help", "Print help messages");
  po::variables_map args;
  try {
    po::store(po::parse_command_line(argc, argv, options), args);
    if (args.count("help")) {
       cerr << options << endl;
       return 0;
    }
    po::notify(args);
  }
  catch (po::error& e) {
    cerr << "ERROR: " << e.what() << endl << endl;
    cerr << options << endl;
    return 1;
  }

  MT19937 eng;
  const string training_corpus_file = args["training_corpus"].as<string>();
  const bool have_initial_alignments = args.count("alignments");
  const string initial_alignment_file = args.count("alignments") ? args["alignments"].as<string>() : "";
  diagonal_alignment_prior diag_alignment_prior(4.0, 0.08, true);
  const unsigned samples = args["samples"].as<int>();
  const bool quiet = args.count("quiet") > 0;
  const bool reverse = args.count("reverse") > 0;
  
  Dict src_dict;
  Dict tgt_dict;
  Dict doc_dict;
  vector<vector<unsigned>> src_corpus;
  vector<vector<unsigned>> tgt_corpus;
  vector<unsigned> document_ids;
  set<unsigned> src_vocab;
  set<unsigned> tgt_vocab;
  ReadFromFile(training_corpus_file, &src_dict, &src_corpus, &src_vocab, &tgt_dict, &tgt_corpus, &tgt_vocab, &doc_dict, &document_ids);
  assert(src_corpus.size() == tgt_corpus.size());
  unsigned document_count = doc_dict.max();
  unsigned tgt_corpus_token_count = 0;
  // dicts contain 1 extra word, <bad>, so the values in src_corpus and tgt_corpus
  // actually run from [1, *_vocab.size()], instead of being 0-indexed.
  cerr << "Corpus size: " << document_count << " documents / " << src_corpus.size() << " sentences\n";
  cerr << src_vocab.size() << " / " << tgt_vocab.size() << " word types\n";

  if (reverse) {
    Dict t1 = src_dict;
    src_dict = tgt_dict;
    tgt_dict = t1;

    vector<vector<unsigned>> t2 = src_corpus;
    src_corpus = tgt_corpus;
    tgt_corpus = t2;

    set<unsigned> t3 = src_vocab;
    src_vocab = tgt_vocab;
    tgt_vocab = t3;
  }

  double uniform_target_word = 1.0 / tgt_vocab.size();
  // Add the null word to the beginning of each source segment
  for (unsigned i = 0; i < src_corpus.size(); ++i) {
    src_corpus[i].insert(src_corpus[i].begin(), 0);
  }

  for (unsigned i = 0; i < tgt_corpus.size(); ++i) {
    tgt_corpus_token_count += tgt_corpus[i].size();
  }

  vector<vector<unsigned short>> alignments;
  vector<vector<set<unsigned short>>> rev_alignments;
  if (have_initial_alignments) {
    assert (!reverse && "Loading alignments with reverse mode not yet implemented");
    cerr << "Loading alignments from " << initial_alignment_file << endl;
    alignments = load_alignment_file(initial_alignment_file, src_corpus, tgt_corpus);
  }

  crp<unsigned> base_ttable(0.0, 1.0);
  vector<crp<unsigned>> underlying_ttable(src_vocab.size() + 1, crp<unsigned>(0.1, uniform_target_word));

  tied_parameter_resampler<crp<unsigned>>        base_ttable_params(1,1,1,1,0.1,1.0);
  tied_parameter_resampler<crp<unsigned>>  underlying_ttable_params(1,1,1,1,0.1,uniform_target_word);

  alignments.resize(tgt_corpus.size());
  rev_alignments.resize(tgt_corpus.size());
  for (unsigned i = 0; i < tgt_corpus.size(); ++i) {
    alignments[i].resize(tgt_corpus[i].size());
    rev_alignments[i].resize(src_corpus[i].size());
  }

  base_ttable_params.insert(&base_ttable);
  for (unsigned i = 0; i < src_vocab.size() + 1; i++) {
    underlying_ttable_params.insert(&underlying_ttable[i]);
  }

  unsigned longest_src_sent_length = 0;
  for (unsigned i = 0; i < src_corpus.size(); i++) {
    longest_src_sent_length = (src_corpus[i].size() > longest_src_sent_length) ? src_corpus[i].size() : longest_src_sent_length;
  } 

  vector<double> a_probs(longest_src_sent_length);
  vector<double> joint_probs(longest_src_sent_length);
  for (unsigned sample=0; sample < samples; ++sample) {
    cerr << "beginning loop with sample = " << sample << endl;
    for (unsigned i = 0; i < tgt_corpus.size(); ++i) {
      const auto& src = src_corpus[i];
      const auto& tgt = tgt_corpus[i];

      for (unsigned n = 0; n < tgt.size(); ++n) {
        unsigned short& a = alignments[i][n];
        const unsigned t = tgt[n];

        if (sample == 0) {
          // random sample during the first iteration
          if (!have_initial_alignments) {
            a = static_cast<unsigned>(sample_uniform01<float>(eng) * src.size());
            if(a == src.size()) {
		cerr << "Warning: sample_uniform01<float>(eng) returned 1.0. Decrementing alignment link by 1." << endl;
                a -= 1;
            }
          }
          assert(a >= 0);
          assert(a < src.size());
        }
        else { 
          rev_alignments[i][a].erase(n);
          if (underlying_ttable[src[a]].decrement(t, eng)) {
            base_ttable.decrement(t, eng);
          }

          // Find the probability of each alignment link

          a_probs.resize(src.size());
          for (unsigned k = 0; k < src.size(); ++k) {
            a_probs[k] = underlying_ttable[src[k]].prob(t, base_ttable.prob(t, uniform_target_word));
            if (k == 0)
              a_probs[k] *= diag_alignment_prior.null_prob(n, tgt.size(), src.size() - 1);
            else
              // TODO: Verify that there's no off by one issue here
              a_probs[k] *= diag_alignment_prior.prob(n + 1, k, tgt.size(), src.size() - 1);
          }

          multinomial_distribution<double> mult(a_probs);
          a = mult(eng);
        }

        // Verify that the draw produced valid results 
        assert(a >= 0);
        assert(a < src.size());
        rev_alignments[i][a].insert(n);

        // Increment the CRPs with the new value
        if (underlying_ttable[src[a]].increment(t, base_ttable.prob(t, uniform_target_word), eng)) {
          base_ttable.increment(t, uniform_target_word, eng);
        }
      }
    }
    output_alignments(tgt_corpus, alignments, reverse);

    if (sample % 10 == 9) {
      cerr << " [LLH=" << log_likelihood(base_ttable_params,
                                         underlying_ttable_params,
                                         diag_alignment_prior,
                                         src_corpus,
                                         tgt_corpus_token_count,
                                         base_ttable,
                                         underlying_ttable,
                                         alignments) << "]" << endl;
      if (sample % 30u == 29) {
        base_ttable_params.resample_hyperparameters(eng);
        underlying_ttable_params.resample_hyperparameters(eng);
        diag_alignment_prior.resample_hyperparameters(alignments, src_corpus, tgt_corpus_token_count, eng);
      }
    }
    else {
      cerr << '.' << flush;
    }
  
    if (sample % 100u == 99 && !quiet) {
      output_latent_variables(underlying_ttable, src_dict, tgt_dict, doc_dict);
    }
  }

  if (true && !quiet) {
    output_latent_variables(underlying_ttable, src_dict, tgt_dict, doc_dict);
  }

  return 0;
}
