#ifndef _OPTIONS_
#define _OPTIONS_
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <getopt.h>
#include "global.h"

/**
 * @brief Prints the help/usage text to CLI
 * 
 */
void print_help_statement();

/**
 * @brief Parses CLI arguments into an Options struct
 * 
 * @param argc arg count
 * @param argv arg vector
 * @param opt output options instance
 */
void parse_options(int argc, char **argv, Options *opt);

#endif /* _OPTIONS_ */
