//
// Created by anonymous on 07.05.2021.
//

#include "graphio.h"

LoadProperties::LoadProperties(const std::map<std::string, int> &propertyMap, const int number) : propertyMap(propertyMap), number(number) {}

bool LoadProperties::checkProperty(std::string path) const {
    for (auto const& [key, value] : propertyMap) {
        if (value == -1){
            if (path.find(key) == std::string::npos) {
                return false;
            }
        }
        else {
            if (path.find(key + std::to_string(value) + "_") == std::string::npos) {
                return false;
            }
        }
    }
    return true;
}
