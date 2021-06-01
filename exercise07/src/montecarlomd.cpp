#include "args.hxx"
#include "mcmd.h"

int main(int argc, char *argv[]) {

  args::ArgumentParser parser("This program performs a molecular dynamics "
                              "simulation using MC in NVE ensemble");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

  args::Group commands(parser, "commands");
  args::Command restart(commands, "restart",
                        "Start from configuration in config.final");
  args::Command equilibration(commands, "equilibration",
                              "Start from configuration in config.0");
  args::Group arguments(parser, "arguments", args::Group::Validators::DontCare,
                        args::Options::Global);
  args::ValueFlag<string> input(arguments, "", "The input file",
                                {'i', "input"});
  args::ValueFlag<string> instantenous(arguments, "",
                                       "yes: print instantaneous values\nno: "
                                       "don't print instantaneous values",
                                       {"instant"});
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    cout << parser;
    return 0;
  } catch (args::ParseError e) {
    cerr << e.what() << endl;
    cerr << parser;
    return 1;
  } catch (args::ValidationError e) {
    cerr << e.what() << endl;
    cerr << parser;
    return 1;
  }

  if (equilibration) {
    if (args::get(instantenous) == "yes") {
      mcmd MCMD(args::get(input), "config.0", true);
      MCMD.Run();
    } else {
      mcmd MCMD(args::get(input), "config.0");
      MCMD.Run();
    }
  }
  if (restart) {
    if (args::get(instantenous) == "yes") {
      mcmd MCMD(args::get(input), "config.final", true);
      MCMD.Run();
    } else {
      mcmd MCMD(args::get(input), "config.final");
      MCMD.Run();
    }
  }
  return 0;
}
