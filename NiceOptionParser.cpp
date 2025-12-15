
#include "NiceOptionParser.hpp"
#include <unistd.h>
#include <stdexcept> 
#include <sstream> 
#include <algorithm> 

using namespace std; 


//______________________________________________________________________________________________________________
NiceOptionParser::NiceOptionParser(int argc, char* argv[], const string& opt_string)
{
    //loop thru the option string. see what options we were handed. 
    int i=0; 
    do {
        //so, we have a new option. let's see what kind it was...
        char opt = opt_string[i]; 

        //these aren't allowed option characters
        if (opt==':' || opt=='?') {
            ostringstream oss; oss << "in <" << __func__ << ">: Invalid option character '" << opt << "'.";
            throw invalid_argument(oss.str()); 
            return;   
        }

        Option_t new_opt{ .flag=opt, .arg_required=false }; 

        //lets see if this option requires an argument (this is the case when the flag is followed by a colon, ex: 'c:')
        if (i < opt_string.size()-1 && opt_string[i+1]==':') {
            i++; 
            new_opt.arg_required=true; 
        }

        //add this option to the list of allowed options
        fOptions_allowed.push_back(new_opt); 

    } while (++i < opt_string.size()); 


    //now, let's actually parse the arguments given. 
    const char* opt_string_char = opt_string.data(); 

    char flag; 
    while ((flag = getopt(argc, argv, opt_string_char)) != -1) {

        //try to find this option in the list of valid options
        auto it = std::find_if( fOptions_allowed.begin(), fOptions_allowed.end(), 
            [flag](const Option_t& option){ return option.flag==flag; }
        ); 

        //this option was NOT found in the list of valid options
        if (flag=='?' || it==fOptions_allowed.end()) {

            ostringstream oss; oss << "in <" << __func__ << ">: Invalid option encountered: '-" << flag << "'"; 
            throw invalid_argument(oss.str()); 
            return; 
        };  

        //we found our option
        Option_t current_option{*it};  

        //check to see if we need an argument for this option
        if (current_option.arg_required) {
            
            if (optarg==nullptr) {
                ostringstream oss; oss << "in <" << __func__ << ">: Option '-" << flag << "' requires argument, but none was provided."; 
                throw invalid_argument(oss.str()); 
                return;
            }

            current_option.arg = string(optarg); 
        }

        fOptions_set.push_back(current_option); 
    }

    //if we've gotten here, then all args were parsed without exception. 
}

//some helper functions...
namespace {

    //returns the iterator correspoding to the option with the given flag
    auto find_option = [](char flag, const vector<NiceOptionParser::Option_t>& opt_vec) {
        auto it = std::find_if(opt_vec.begin(), opt_vec.end(), [flag](const NiceOptionParser::Option_t& opt){ return opt.flag==flag; }); 
        return it; 
    };

    //returns 'true' if a flag exists in the vector of possible options, and 'false' otherwise
    auto option_exists = [](char flag, const vector<NiceOptionParser::Option_t>& opt_vec) {
        auto it = find_option(flag, opt_vec); 
        return (it != opt_vec.end()); 
    };
}

//______________________________________________________________________________________________________________
bool NiceOptionParser::Is_option_set(char flag) const 
{
    //see if this option is allowed
    if (!option_exists(flag, fOptions_allowed)) {

        ostringstream oss; oss << "in <" << __func__ << ">: Invalid option: '-" << flag << "'"; 
        throw invalid_argument(oss.str()); 
    }

    //now that we know that 'flag' is an allowed option (that may or may not be set), 
    // we return 'true' if it was invoked on this execution, and 'false' if it was not. 
    return option_exists(flag, fOptions_set); 
}

//______________________________________________________________________________________________________________
string NiceOptionParser::Arg_or_defalut_value(char flag, string default_val) const 
{
    //see if this option is allowed
    if (!option_exists(flag, fOptions_allowed)) {

        ostringstream oss; oss << "in <" << __func__ << ">: Invalid option: '-" << flag << "'"; 
        throw invalid_argument(oss.str()); 
    }

    Option_t opt = *find_option(flag, fOptions_allowed); 
    
    //see if this flag actually requires an argument
    if (!opt.arg_required) {
        
        ostringstream oss; oss << "in <" << __func__ << ">: Option '-" << flag << "' does not take arguments"; 
        throw invalid_argument(oss.str()); 
    }
    
    //now, let's see if we should return the argument or the default value
    if (option_exists(flag, fOptions_set)) {

        return find_option(flag, fOptions_set)->arg; 
    
    } else {

        return default_val; 
    }
}