#ifndef NiceOptionParser_HPP
#define NiceOptionParser_HPP

#include <string> 
#include <vector>  

class NiceOptionParser {
public: 
    struct Option_t { 
        char flag;
        bool arg_required=false; 
        std::string arg{""}; 
        
        //returns true if these are the same options
        bool operator==(const Option_t& rhs) const { return flag==rhs.flag; };  
    }; 

    //Create the NiceOptionParser object, and feed it the list of options/args
    NiceOptionParser(int argc, char* argv[], const std::string& opt_string); 

    //has a particular option been set? 
    bool Is_option_set(char opt) const; 

    //in the style of 'std::optional::value_or()', this returns either the argument passed to a particular option, 
    // or it returns the 'default_val' value provided if the user did not specify this flag on the invocation of 'main'. 
    std::string Arg_or_defalut_value(char opt, std::string default_val) const; 

private: 
    //set of all allowed options
    std::vector<Option_t> fOptions_allowed; 

    //set of all options that were invoked on this particular execution
    std::vector<Option_t> fOptions_set;

}; 


#endif 