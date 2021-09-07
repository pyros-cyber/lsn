#include "TSP.h"
#include "args.hxx"

int main(int argc, char *argv[]) {
  args::ArgumentParser parser(
      "This program solves the TSP usign a generic algorithm");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Group arguments(parser, "arguments", args::Group::Validators::All,
                        args::Options::Global);
  args::ValueFlag<double> Npop(
      arguments, "", "Number of individuals in the first generation", {"Npop"});
  args::ValueFlag<double> Ngen(arguments, "", "Number of generations",
                               {"Ngen"});
  args::ValueFlag<string> shape(
      arguments, "",
      "Shape of the world, can be either circumference or square", {"shape"});
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
  // double cross = 0.7;
  Random rnd("Primes", "seed.in");
  int _Ncities = 32;
  int _Npop = args::get(Npop);
  int _Ngen = args::get(Ngen);
  string _shape = args::get(shape);

  ofstream cities("./results/cities_" + _shape + ".dat");
  ofstream res("./results/results_" + _shape + ".txt");
  if (cities.is_open() && res.is_open()) {
    Population pop(make_shared<Random>(rnd), _shape, _Npop, _Ncities);

    for (int i{}; i < _Ncities; ++i) {
      cities << pop.Pop[0].Cities[pop.Pop[0].Path[i] - 1].first << " "
             << pop.Pop[0].Cities[pop.Pop[0].Path[i] - 1].second << endl;
    }
    cities.close();

    cities.open("./results/cities_final_" + _shape + ".dat");
    if (cities.is_open()) {
      pop.OrderPop();

      cout << "------------------------------" << endl;
      cout << _Ncities << " cities on a " << _shape << ";" << endl;
      cout << "Population of " << _Npop << " salesmen;" << endl;
      cout << "Evolution of " << _Ngen << " generations. \n" << endl;
      cout << "Best Initial loss value: " << pop.Losses[_Npop - 1] << endl;

      res << "------------------------------" << endl;
      res << _Ncities << " cities on a " << _shape << ";" << endl;
      res << "Population of " << _Npop << " salesmen;" << endl;
      res << "Evolution of " << _Ngen << " generations. \n" << endl;
      res << "Best Initial loss value: " << pop.Losses[_Npop - 1] << endl;

      pop.Evolve(_Ngen, _shape);

      cout << "Best final loss value: " << pop.Losses[_Npop - 1] << endl;
      res << "Best final loss value: " << pop.Losses[_Npop - 1] << endl;

      for (int i{}; i < _Ncities; ++i) {
        cities
            << pop.Pop[_Npop - 1].Cities[pop.Pop[_Npop - 1].Path[i] - 1].first
            << "  "
            << pop.Pop[_Npop - 1].Cities[pop.Pop[_Npop - 1].Path[i] - 1].second
            << endl;
      }

      cities << pop.Pop[_Npop - 1].Cities[pop.Pop[_Npop - 1].Path[0] - 1].first
             << "  "
             << pop.Pop[_Npop - 1].Cities[pop.Pop[_Npop - 1].Path[0] - 1].second
             << endl;

      cities.close();
      //  	losses.close();
      //  	losses_ave.close();

    } else {
      cerr << "ERROR: can't open output files. Exiting." << endl;
      return 1;
    }
  } else {
    cerr << "ERROR: can't open output files. Exiting." << endl;
    return 1;
  }
  return 0;
}
