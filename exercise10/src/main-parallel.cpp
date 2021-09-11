#include "TSP.h"
#include "args.hxx"

int main(int argc, char *argv[]) {
  args::ArgumentParser parser("This program solves the TSP using a parall "
                              "genetic algorithm");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Group arguments(parser, "arguments", args::Group::Validators::All,
                        args::Options::Global);
  args::ValueFlag<int> Npop(arguments, "",
                            "GA: Number of individuals in the first generation",
                            {"Npop"});
  args::ValueFlag<int> Ngen(arguments, "", "GA: Number of generations",
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
  int Rank;
  int size;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  MPI_Status stat;

  Random rnd(Rank);

  /* we are using it in std::random_shuffle,
   * it has to be different for each run of the program!
   */
  srand(Rank + 1);

  int _Ncities = 32;
  int _Nmigr = 10;
  int _Npop = args::get(Npop);
  int _Ngen = args::get(Ngen);
  string _shape = args::get(shape);

  ofstream cities("./results/cities_initial_" + _shape + "_" + to_string(Rank) +
                  ".dat");
  if (cities.is_open()) {
    Population pop(make_shared<Random>(rnd), _shape, _Npop, _Ncities);
    vector<pair<double, double>> _Cities(_Ncities);
    vector<double> x_coord(_Ncities);
    vector<double> y_coord(_Ncities);

    if (Rank == 0) {
      for (int i{}; i < _Ncities; ++i) {
        x_coord[i] = pop.Pop[0].GetCities()[i].first;
        y_coord[i] = pop.Pop[0].GetCities()[i].second;
      }
    }

    MPI_Bcast(&x_coord.front(), _Ncities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&y_coord.front(), _Ncities, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i{}; i < _Ncities; ++i) {
      _Cities[i].first = x_coord[i];
      _Cities[i].second = y_coord[i];
    }

    for (int i{}; i < _Npop; ++i) {
      pop.Pop[i].Cities = _Cities;
    }

    for (int i{}; i < _Ncities; ++i) {
      cities << pop.Pop[0].GetCities()[i].first << " "
             << pop.Pop[0].GetCities()[i].second << endl;
    }

    cities << pop.Pop[0].GetCities()[0].first << "  "
           << pop.Pop[0].GetCities()[0].second << endl;

    cities.close();

    cities.open("./results/cities_final_" + _shape + "_" + to_string(Rank) +
                ".dat");

    if (cities.is_open()) {
      pop.OrderPop();

      if (Rank == 0) {
        cout << "------------------------------" << endl;
        cout << _Ncities << " cities on a " << _shape << ";" << endl;
        cout << "Population of " << _Npop << " salesmen;" << endl;
        cout << "Evolution of " << _Ngen << " generations. \n" << endl;
      }

      cout << "Rank_" << Rank
           << ") Best Initial loss value: " << pop.Losses[_Npop - 1] << endl;

      pop.ParallelEvolve(Rank, _Ngen, _Nmigr, _shape, &stat);

      cout << "Rank_" << Rank
           << ") Best final loss value: " << pop.Losses[_Npop - 1] << endl;

      for (int i{}; i < _Ncities; ++i) {
        cities << pop.Pop[_Npop - 1].GetCities()[i].first << "  "
               << pop.Pop[_Npop - 1].GetCities()[i].second << endl;
      }

      cities << pop.Pop[_Npop - 1].GetCities()[0].first << "  "
             << pop.Pop[_Npop - 1].GetCities()[0].second << endl;

      cities.close();
      MPI_Finalize();
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
