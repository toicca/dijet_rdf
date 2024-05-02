#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>
#include <string>

using json = nlohmann::json;
    
json golden_json;

void init_json(std::string jsonFile) {
    std::cout << "Initializing JSON file" << std::endl;
    std::ifstream f(jsonFile);
    golden_json = json::parse(f);
}

bool isGoodLumi(int run, int lumi) {    
   for (auto& lumiRange : golden_json[std::to_string(run)]) {
       if (lumi >= lumiRange[0] && lumi <= lumiRange[1]) {
           return true;
       }
   }
    
    return false;
}