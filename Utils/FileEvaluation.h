//
// Created by anonymous on 21.10.21.
//

#ifndef CLOSURES_FILEEVALUATION_H
#define CLOSURES_FILEEVALUATION_H


#include <string>
#include <vector>
#include <iomanip>
#include <unordered_map>

class FileEvaluation {
    std::string out_path;
    std::string name;
    std::vector<std::string> headers_all;
    std::unordered_map<std::string, std::vector<std::string>> values_all;
    std::vector<std::string> headers_summary;
    std::unordered_map<std::string, std::vector<std::string>> values_summary;
    std::string extension = ".csv";
public:
    FileEvaluation(){};
    explicit FileEvaluation(const std::string& out_path, const std::string& eval_name = "default", const std::string& extension = ".csv") : out_path(out_path), name(eval_name), extension(extension){};
    void save(bool summary = false, bool both = true, std::_Ios_Openmode mode = std::ios_base::app);
    void headerValueInsert(const std::vector<std::string>& new_header, const std::vector<std::string>& new_values, int insert_position=-1, bool summary = false, bool both=false);
    void clear();

};


#endif //CLOSURES_FILEEVALUATION_H
