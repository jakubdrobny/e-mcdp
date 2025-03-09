#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include <string>

class Output {
public:
  Output(const std::string &filename);
  ~Output();
  void print(const std::string &msg);

private:
  std::ofstream out;
  bool out_to_file;
  std::string filename;
};

#endif // OUTPUT_H
