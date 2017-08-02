#include <iostream>
#include <argp.h>

#include "InputFileReader.h"
#include "Config.h"

using std::cout;
using std::endl;

//const char *sean_program_version = "SeAN 0.0.0";
//const char *sean_program_bug_address = "<ugayer@ikp.tu-darmstadt.de>";

static char doc[] = "SeAN, Self-Absorption Numerical";
static char args_doc[] = "INPUTFILE";

static struct argp_option options[] = {
  { 0, 0, 0, 0, 0 }
};

struct arguments{
        char *inputfile;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = (struct arguments*)state->input;

    switch (key) {
    case ARGP_KEY_ARG: arguments->inputfile = arg; break;
    case ARGP_KEY_END:
        if(state->arg_num == 0) {
            argp_usage(state);
        }
        if(!arguments->inputfile) {
            argp_failure(state, 1, 0, "INPUTFILE not specified.");
        }
        break;
    default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

int main(int argc, char* argv[]){

	struct arguments arguments;
	argp_parse(&argp, argc, argv, 0, 0, &arguments);

	InputFileReader inputFileReader;

	inputFileReader.readInputFile(arguments.inputfile);
	inputFileReader.print();
}
