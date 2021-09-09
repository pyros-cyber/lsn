#include "TSP.h"
#include "args.hxx"

int main(int argc, char *argv[]) {
  args::ArgumentParser parser("This program solves the TSP using both a "
                              "genetic and simulated annealing algorithms");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Group arguments(parser, "arguments", args::Group::Validators::All,
                        args::Options::Global);
  args::ValueFlag<int> Npop(arguments, "",
                            "GA: Number of individuals in the first generation",
                            {"Npop"});
  args::ValueFlag<int> Ngen(arguments, "", "GA: Number of generations",
                            {"Ngen"});
  args::ValueFlag<int> ntemp(arguments, "", "SA: Temperatures to consider",
                             {"ntemp"});
  args::ValueFlag<int> nsteps(arguments, "", "SA: MC steps of the simulation",
                              {"nsteps"});
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
  int _ntemp = args::get(ntemp);
  int _nsteps = args::get(nsteps);
  string _shape = args::get(shape);
  vector<pair<double, double>> initial_coordinates;

  ofstream cities("./results/ga_cities_" + _shape + ".dat");
  ofstream res("./results/ga_results_" + _shape + ".txt");
  if (cities.is_open() && res.is_open()) {
    Population pop(make_shared<Random>(rnd), _shape, _Npop, _Ncities);
    initial_coordinates = pop.Pop[0].GetCities();

    for (int i{}; i < _Ncities; ++i) {
      cities << pop.Pop[0].GetCities()[i].first << " "
             << pop.Pop[0].GetCities()[i].second << endl;
    }
    cities.close();

    cities.open("./results/ga_cities_final_" + _shape + ".dat");
    if (cities.is_open()) {
      pop.OrderPop();

      cout << "------------------------------" << endl;
      cout << "Genetic Algorithm" << endl;
      cout << _Ncities << " cities on a " << _shape << ";" << endl;
      cout << "Population of " << _Npop << " salesmen;" << endl;
      cout << "Evolution of " << _Ngen << " generations. \n" << endl;
      cout << "Best initial loss value: " << pop.Losses[_Npop - 1] << endl;

      res << "------------------------------" << endl;
      res << "Genetic Algorithm" << endl;
      res << _Ncities << " cities on a " << _shape << ";" << endl;
      res << "Population of " << _Npop << " salesmen;" << endl;
      res << "Evolution of " << _Ngen << " generations. \n" << endl;
      res << "Best initial loss value: " << pop.Losses[_Npop - 1] << endl;

      pop.Evolve(_Ngen, _shape);

      cout << "Best final loss value: " << pop.Losses[_Npop - 1] << endl;
      res << "Best final loss value: " << pop.Losses[_Npop - 1] << endl;

      for (int i{}; i < _Ncities; ++i) {
        cities << pop.Pop[_Npop - 1].GetCities()[i].first << "  "
               << pop.Pop[_Npop - 1].GetCities()[i].second << endl;
      }

      cities << pop.Pop[_Npop - 1].GetCities()[0].first << "  "
             << pop.Pop[_Npop - 1].GetCities()[0].second << endl;

      cities.close();
      res.close();

    } else {
      cerr << "ERROR: can't open output files. Exiting." << endl;
      return 1;
    }
  } else {
    cerr << "ERROR: can't open output files. Exiting." << endl;
    return 1;
  }

  cities.open("./results/sa_cities_" + _shape + ".dat");
  res.open("./results/sa_results_" + _shape + ".txt");

  if (cities.is_open() && res.is_open()) {
    Salesman Initial(make_shared<Random>(rnd), _shape, _Ncities);
    Initial.Cities = initial_coordinates;
    Salesman Final = Initial;

    for (int i{}; i < _Ncities; ++i) {
      cities << Initial.GetCities()[i].first << " "
             << Initial.GetCities()[i].second << endl;
    }

    cities.close();
    cities.open("./results/sa_cities_final_" + _shape + ".dat");

    if (cities.is_open()) {
      cout << "------------------------------" << endl;
      cout << "Simulated Annealing Algorithm" << endl;
      cout << _Ncities << " cities on a " << _shape << ";" << endl;
      cout << _ntemp << " temperatures to be considered;" << endl;
      cout << _nsteps << " steps to be made for each temperature. \n" << endl;
      cout << "Best initial loss value: " << Initial.AbsoluteLoss() << endl;

      res << "------------------------------" << endl;
      res << "Simulated Annealing Algorithm" << endl;
      res << _Ncities << " cities on a " << _shape << ";" << endl;
      res << _ntemp << " temperatures to be considered;" << endl;
      res << _nsteps << " steps to be made for each temperature. \n" << endl;
      res << "Best initial loss value: " << Initial.AbsoluteLoss() << endl;

      Final.SimulatedAnnealing(Initial, _shape, _ntemp, _nsteps);

      cout << "Final loss value: " << Initial.AbsoluteLoss() << endl;
      res << "Final loss value: " << Initial.AbsoluteLoss() << endl;

      for (int i{}; i < _Ncities; ++i) {
        cities << Initial.GetCities()[i].first << " "
               << Initial.GetCities()[i].second << endl;
      }

      cities << Initial.GetCities()[0].first << " "
             << Initial.GetCities()[0].second << endl;

      cities.close();
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
