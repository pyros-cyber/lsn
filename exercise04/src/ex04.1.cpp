#include "MolDyn_NVE.h"
#include "args.hxx"

using namespace std;

int interface(int argc, char *argv[]) {
  args::ArgumentParser parser(
      "This program performs a molecular dynamics simulation in NVE ensemble");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<string> restart(
      parser, "", "Start from configuration in old.0 and config.0",
      {'r', "restart"});
  args::ValueFlag<string> equilibrium(
      parser, "", "Start from configuration in config.0", {'e', "equilibrium"});
  args::ValueFlag<string> input(parser, "", "The input file", {'i', "input"});

  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    cout << parser;
    return 0;
  } catch (args::ParseError e) {
    cerr << e.what() << std::endl;
    cerr << parser;
    return 1;
  } catch (args::ValidationError e) {
    cerr << e.what() << std::endl;
    cerr << parser;
    return 1;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  interface(argc, argv);
  return 0;
}