#include "args.hxx"
#include "hydrogen.h"
#include "metropolis.h"
#include "myStatFunc.h"
#include "random.h"

using namespace std;

int main(int argc, char *argv[]) {

  args::ArgumentParser parser(
      "This program estimates the value of <r> for the ground state and first "
      "excited state of the H2 atom");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<string> state(
      parser, "", "State of H2 atom, can be either 'ground' or 'first-excited'",
      {'s', "state"});
  args::ValueFlag<string> transition(
      parser, "",
      "Transition probability to use in metropolis algorithm. Can either be "
      "'uniform' or 'gaussian'",
      {'t', "transition"});
  args::ValueFlag<double> step(parser, "", "Step of the Metropolis algorithm",
                               {"step"});
  args::ValueFlag<double> x_position(
      parser, "", "x coordinate of the starting position", {'x', "x_pos"});
  args::ValueFlag<double> y_position(
      parser, "", "y coordinate of the starting position", {'y', "y_pos"});
  args::ValueFlag<double> z_position(
      parser, "", "z coordinate of the starting position", {'z', "z_pos"});
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

  unsigned int M = 1e7, N = 100;
  int throws_per_block = M / N;

  HydrogenAtom hydro(args::get(x_position), args::get(y_position),
                     args::get(z_position), args::get(state));
  Metropolis metro(args::get(step), hydro.x, hydro.y, hydro.z);
  vector<double> r(N, 0.);
  for (int i{}; i < N; ++i) {
    for (int j{}; j < throws_per_block; ++j) {

      if (args::get(transition) == "uniform" && hydro.state == "ground") {
        metro.Run(
            [&metro](double coordinate) {
              return metro.uniform_step(coordinate);
            },
            ground_state);
      } else if (args::get(transition) == "uniform" &&
                 hydro.state == "first-excited") {
        metro.Run(
            [&metro](double coordinate) {
              return metro.uniform_step(coordinate);
            },
            first_excited_state);
      } else if (args::get(transition) == "gaussian" &&
                 hydro.state == "ground") {
        metro.Run(
            [&metro](double coordinate) {
              return metro.gaussian_step(coordinate);
            },
            ground_state);
      } else if (args::get(transition) == "gaussian" &&
                 hydro.state == "first-excited") {
        metro.Run(
            [&metro](double coordinate) {
              return metro.gaussian_step(coordinate);
            },
            first_excited_state);
      } else {
        cerr << "Unknown transition probability or hydrogen state.\nExiting.\n";
        return 1;
      }

      r[i] +=
          sqrt(metro.x_start * metro.x_start + metro.y_start * metro.y_start +
               metro.z_start * metro.z_start);
    }
    r[i] /= throws_per_block;
  }

  cout << "Metropolis algorithm acceptance: " << metro.acceptance() << endl;

  vector<double> radius_error = blocking_error(r);

  ofstream out("results." + args::get(transition) + "." + hydro.state + ".dat");
  if (out.is_open()) {
    for (int i{}; i < r.size(); ++i) {
      out << r[i] << " " << radius_error[i] << endl;
    }
  } else {
    cerr << "ERROR: unable to open output file" << endl;
    return 1;
  }
  out.close();
  return 0;
}