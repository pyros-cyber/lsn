#include "MolDyn_NVE.h"
#include "args.hxx"

using namespace std;

int interface(int argc, char *argv[]) {
  args::ArgumentParser parser(
      "This program performs a molecular dynamics simulation in NVE ensemble");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

  args::Group commands(parser, "commands");
  args::Command restart(commands, "restart",
                        "Start from configuration in old.0 and config.final");
  args::Command equilibration(commands, "equilibration",
                              "Start from configuration in config.0");
  args::Group arguments(parser, "arguments", args::Group::Validators::DontCare,
                        args::Options::Global);
  args::ValueFlag<string> input(arguments, "", "The input file",
                                {'i', "input"});

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
  if (restart) {
    MolDyn_NVE molecularDynamic(args::get(input), "config.final", "old.0");
    molecularDynamic.RunSimulation();
  }
  if (equilibration) {
    MolDyn_NVE molecularDynamic(args::get(input), "config.0");
    molecularDynamic.RunSimulation();
  }
  return 0;
}

int main(int argc, char *argv[]) {
  interface(argc, argv);
  return 0;
}