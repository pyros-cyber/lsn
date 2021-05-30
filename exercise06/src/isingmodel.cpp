#include "args.hxx"
#include "ising.h"
#include "random.h"

using namespace std;

int main(int argc, char *argv[]) {

  args::ArgumentParser parser(
      "This program simulates the behaviour of a 1D Ising model, using either "
      "a metropolis algorithm or a gibbs algorithm");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Group commands(parser, "commands");
  args::Command equilibration(
      commands, "equilibration",
      "Start from random configuration and print "
      "instantaneous values of the observed quantities");
  args::Command first(commands, "first", "Start from random configuration");
  args::Command restart(commands, "restart",
                        "Restart configuration from file 'config.final'");
  args::Group arguments(parser, "arguments", args::Group::Validators::DontCare,
                        args::Options::Global);
  args::ValueFlag<string> algorithm(
      arguments, "", "Algorithm to use, can either be 'metropolis' or 'gibbs'",
      {'a', "algorithm"});
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
  if (equilibration) {
    Ising1D isingmodel;
    if (args::get(algorithm) == "metropolis") {
      isingmodel.Run([&isingmodel]() { return isingmodel.MetropolisMove(); },
                     true);
    } else if (args::get(algorithm) == "gibbs") {
      isingmodel.Run([&isingmodel]() { return isingmodel.GibbsMove(); }, true);
    } else {
      cerr << "Unknown algorithm.\nExiting.\n";
      return 1;
    }
  }
  if (first) {
    Ising1D isingmodel;
    if (args::get(algorithm) == "metropolis") {
      isingmodel.Run([&isingmodel]() { return isingmodel.MetropolisMove(); });
    } else if (args::get(algorithm) == "gibbs") {
      isingmodel.Run([&isingmodel]() { return isingmodel.GibbsMove(); });
    } else {
      cerr << "Unknown algorithm.\nExiting.\n";
      return 1;
    }
  }
  if (restart) {
    Ising1D isingmodel("config.final");
    if (args::get(algorithm) == "metropolis") {
      isingmodel.Run([&isingmodel]() { return isingmodel.MetropolisMove(); });
    } else if (args::get(algorithm) == "gibbs") {
      isingmodel.Run([&isingmodel]() { return isingmodel.GibbsMove(); });
    } else {
      cerr << "Unknown algorithm.\nExiting.\n";
      return 1;
    }
  }

  return 0;
}