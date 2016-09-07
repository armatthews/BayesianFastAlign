#ifndef CPYPDICT_H_
#define CPYPDICT_H_

#include <string>
#include <iostream>
#include <cassert>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_map>
#include <functional>

namespace cpyp {

class Dict {
 typedef std::unordered_map<std::string, unsigned, std::hash<std::string> > Map;
 public:
  Dict() : b0_("<bad0>") {
    words_.reserve(1000);
  }

  inline unsigned max() const { return words_.size(); }

  static bool is_ws(char x) {
    return (x == ' ' || x == '\t');
  }

  inline void ConvertWhitespaceDelimitedLine(const std::string& line, std::vector<unsigned>* out) {
    size_t cur = 0;
    size_t last = 0;
    int state = 0;
    out->clear();
    while(cur < line.size()) {
      if (is_ws(line[cur++])) {
        if (state == 0) continue;
        out->push_back(Convert(line.substr(last, cur - last - 1)));
        state = 0;
      } else {
        if (state == 1) continue;
        last = cur - 1;
        state = 1;
      }
    }
    if (state == 1)
      out->push_back(Convert(line.substr(last, cur - last)));
  }

  inline unsigned Convert(const std::string& word, bool frozen = false) {
    Map::iterator i = d_.find(word);
    if (i == d_.end()) {
      if (frozen)
        return 0;
      words_.push_back(word);
      d_[word] = words_.size();
      return words_.size();
    } else {
      return i->second;
    }
  }

  inline const std::string& Convert(const unsigned id) const {
    if (id == 0) return b0_;
    return words_[id-1];
  }
  template<class Archive> void serialize(Archive& ar, const unsigned int version) {
    ar & b0_;
    ar & words_;
    ar & d_;
  }
 private:
  std::string b0_;
  std::vector<std::string> words_;
  Map d_;
};

std::vector<std::string> split(const std::string& s) {
  std::vector<std::string> result;
  unsigned j = 0;
  for(unsigned i = 0; i < s.length() - 2; i++) {
    if(s[i] == '|' && s[i+1] == '|' && s[i+2] == '|') {
      if(i != j) {
        result.push_back(s.substr(j, i - j));
        j = i + 3;
      }
    }
  }
  if(j != s.length()) {
    result.push_back(s.substr(j));
  }
  return result;
}

void ReadFromFile(const std::string& filename,
                  Dict* src_dict,
                  std::vector<std::vector<unsigned> >* src,
                  std::set<unsigned>* src_vocab,
                  Dict* tgt_dict,
                  std::vector<std::vector<unsigned> >* tgt,
                  std::set<unsigned>* tgt_vocab,
                  Dict* doc_dict = NULL,
                  std::vector<unsigned>* doc_ids = NULL) {
  src->clear();
  std::cerr << "Reading from " << filename << std::endl;
  std::ifstream in(filename);
  assert(in);
  std::string line;
  int lc = 0;
  while(getline(in, line)) {
    ++lc;
    std::vector<std::string> parts = split(line);

    std::string src_line, tgt_line, doc_id;
    if(parts.size() == 2) {
      doc_id = "";
      src_line = parts[0];
      tgt_line = parts[1];
    }
    else {
      doc_id = parts[0];
      src_line = parts[1];
      tgt_line = parts[2];
    }

    src->push_back(std::vector<unsigned>());
    src_dict->ConvertWhitespaceDelimitedLine(src_line, &src->back());
    for (unsigned i = 0; i < src->back().size(); ++i) src_vocab->insert(src->back()[i]);

    tgt->push_back(std::vector<unsigned>());
    tgt_dict->ConvertWhitespaceDelimitedLine(tgt_line, &tgt->back());
    for (unsigned i = 0; i < tgt->back().size(); ++i) tgt_vocab->insert(tgt->back()[i]);

    if (doc_ids != NULL && doc_dict != NULL)
      doc_ids->push_back(doc_dict->Convert(doc_id));
  }
}

}

#endif
