#include <iostream>
#include <sstream>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>
using namespace std;
namespace po = boost::program_options;

vector<string> tokenize(string input, string delimiter, unsigned max_times);
vector<string> tokenize(string input, string delimiter);
vector<string> tokenize(string input, char delimiter);

string strip(const string& input);
vector<string> strip(const vector<string>& input, bool removeEmpty = false);

vector<string> tokenize(string input, string delimiter, unsigned max_times) {
  vector<string> tokens;
  //tokens.reserve(max_times);
  size_t last = 0;
  size_t next = 0;
  while ((next = input.find(delimiter, last)) != string::npos && tokens.size() < max_times) {
    tokens.push_back(input.substr(last, next-last));
    last = next + delimiter.length();
  }
  tokens.push_back(input.substr(last));
  return tokens;
}

vector<string> tokenize(string input, string delimiter) {
  return tokenize(input, delimiter, input.length());
}

vector<string> tokenize(string input, char delimiter) {
  return tokenize(input, string(1, delimiter));
}

string strip(const string& input) {
  string output = input;
  boost::algorithm::trim(output);
  return output;
}

vector<string> strip(const vector<string>& input, bool removeEmpty) {
  vector<string> output;
  for (unsigned i = 0; i < input.size(); ++i) {
    string s = strip(input[i]);
    if (s.length() > 0 || !removeEmpty) {
      output.push_back(s);
    }
  }
  return output;
}

vector<pair<unsigned, unsigned>> ParseLine(const string& line) {
  vector<pair<unsigned, unsigned>> links;
  vector<string> link_strings = tokenize(line, ' ');
  for (const string& link_string : link_strings) {
    if (link_string.length() == 0) {
      continue;
    }
    vector<string> parts = tokenize(link_string, '-');
    assert (parts.size() == 2);
    unsigned i = atoi(parts[0].c_str());
    unsigned j = atoi(parts[1].c_str());
    links.push_back(make_pair(i, j));
  }
  return links;
}

int main(int argc, char** argv) {
  po::options_description desc("description");
  desc.add_options()
  ("help", "Display this help message")
  ("burn_in,b", po::value<unsigned>()->default_value(0), "Training text, morphologically analyzed")
  ("reverse,r", "Compute mode w.r.t. each source word rather than each target word. Use this iff you used -r with bayesianfastalign.")
  ("sentence_count,n", po::value<unsigned>()->required(), "Number of sentences in the training corpus");

  po::positional_options_description positional_options;
  positional_options.add("sentence_count", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(positional_options).run(), vm);

  if (vm.count("help")) {
    cerr << desc;
    return 1;
  }

  po::notify(vm);

  const unsigned sentence_count = vm["sentence_count"].as<unsigned>();
  const unsigned burn_in = vm["burn_in"].as<unsigned>();
  const bool reverse = vm.count("reverse") > 0;

  vector<vector<vector<unsigned>>> link_counts(sentence_count); // link_counts[sent_index][target_index][source_index] = count

  unsigned line_num = 0;
  unsigned iteration = 0;
  ios_base::sync_with_stdio (false); // supposedly this makes reading from stdin faster
  for (string line; getline(cin, line); ++line_num) {
    unsigned sent_num = line_num % sentence_count;
    if (sent_num == 0 && line_num != 0) {
      iteration++;
    }

    if (iteration < burn_in) {
      continue;
    }

    char* p = (char*)line.c_str();
    unsigned i, j;
    int n;
    while (sscanf(p, "%u-%u%n", &i, &j, &n) != -1) {
      if (reverse) {
        unsigned t = i;
        i = j;
        j = t;
      }
      assert (p < line.c_str() + line.length());
      p += n;

      assert (sent_num < link_counts.size());
      if (j >= link_counts[sent_num].size()) {
        link_counts[sent_num].resize(j + 1);
      }

      assert (j < link_counts[sent_num].size());
      if (i >= link_counts[sent_num][j].size()) {
        link_counts[sent_num][j].resize(i + 1, 0);
      }

      assert (i < link_counts[sent_num][j].size());
      link_counts[sent_num][j][i]++;
    }
  }

  for (unsigned n = 0; n < sentence_count; ++n) {
    bool first = true;
    for (unsigned j = 0; j < link_counts[n].size(); ++j) {
      if (link_counts[n][j].size() == 0) {
        continue;
      }
      unsigned best_i = 0;
      unsigned best_v = link_counts[n][j][0];
      unsigned sum = best_v;
      for (unsigned i = 1; i < link_counts[n][j].size(); ++i) {
        unsigned v = link_counts[n][j][i];
        sum += v;
        if (v > best_v) {
          best_v = v;
          best_i = i;
        }
      }

      if (best_v >= (iteration + 1 - burn_in) - sum) {
        if (!first) {
          cout << " ";
        }
        else {
          first = false;
        }
        if (reverse) {
          cout << j << "-" << best_i;
        }
        else {
          cout << best_i << "-" << j;
        }
      }
    }
    cout << endl;
  }

  return 0;
}
