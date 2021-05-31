#include "mcmd.h"

mcmd::mcmd(string simParameters, string initial_configuration,
           bool _instantenous)
    : rnd{"../Primes", "../seed.in"}, nbins{100} {
  instantenous = _instantenous;
  props = {"energy", "pressure"};
  for (auto &el : props) {
    walker[el] = 0.;
    block_average[el] = 0.;
    glob_average[el] = 0.;
    glob_average_sq[el] = 0.;
  }
  g_histo.resize(nbins);
  g_histo_block_ave.resize(nbins);
  g_histo_glob_ave.resize(nbins);
  g_histo_glob_ave_sq.resize(nbins);

  blk_norm = 0;
  accepted = 0;
  attempted = 0;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl
       << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T "
       << endl
       << endl;
  cout << "The program uses Lennard-Jones units \n" << endl;
  cout << "Reading parameters from " << simParameters << endl;

  /* READING INPUT FILES AND INITIALIZING PARAMETERS */
  ifstream ReadInput(simParameters);
  if (ReadInput.fail()) {
    cerr << "ERROR: unable to open " << simParameters << " file.\nExiting.\n";
    exit(1);
  }
  ReadInput >> temp;
  beta = 1.0 / temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = static_cast<double>(npart) / rho;
  box = pow(vol, 1. / 3.);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;

  // Tail corrections for potential energy and pressure
  vtail = (8.0 * M_PI * rho) / (9.0 * pow(rcut, 9)) -
          (8.0 * M_PI * rho) / (3.0 * pow(rcut, 3));
  ptail = (32.0 * M_PI * rho) / (9.0 * pow(rcut, 9)) -
          (16.0 * M_PI * rho) / (3.0 * pow(rcut, 3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the pressure         = " << ptail << endl;

  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblk;

  cout << "The program perform Metropolis moves with uniform translations"
       << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

  bin_size = (box / 2.) / nbins;

  x.resize(npart);
  y.resize(npart);
  z.resize(npart);

  ifstream ReadConf(initial_configuration);
  if (ReadConf.fail()) {
    cerr << "ERROR: unable to open " << initial_configuration
         << " file.\nExiting.\n";
    exit(1);
  }
  for (int i{}; i < npart; ++i) {
    double x_temp, y_temp, z_temp;
    ReadConf >> x_temp >> y_temp >> z_temp;
    x[i] = Pbc(x_temp * box);
    y[i] = Pbc(y_temp * box);
    z[i] = Pbc(z_temp * box);
  }
  ReadConf.close();

  // Evaluate potential energy and pressure of the initial configuration
  Measure();
  // Print initial values for the potential energy and pressure
  cout << "Initial potential energy (with tail corrections) = "
       << walker.at("energy") / npart + vtail << endl;
  cout << "Virial                   (with tail corrections) = "
       << walker.at("pressure") / npart + ptail << endl;
  cout << "Pressure                 (with tail corrections) = "
       << rho * temp + (walker.at("pressure") + npart * ptail) / vol << endl
       << endl;

  /* OPENING OUTPUT STREAMS */
  Epot.open("results/epot.dat");
  Pres.open("results/pres.dat");
  Gerr.open("results/gerr.dat");
  Gave.open("results/gave.dat");
  Binn.open("results/binning.dat");

  if (instantenous) {
    ist_pot.open("results/instant_epot.dat");
    ist_pres.open("results/instant_pres.dat");
  }
}

mcmd::~mcmd() {
  Epot.close();
  Pres.close();
  Gerr.close();
  Gave.close();
  if (instantenous) {
    ist_pot.close();
    ist_pres.close();
  }
}

double mcmd::Boltzmann(double _x, double _y, double _z, unsigned int ip) {
  double ene = 0., dx = 0., dy = 0., dz = 0., dr = 0.;

  for (int i{}; i < npart; ++i) {
    if (i != ip) {
      // distance ip-i in pbc
      dx = Pbc(_x - x[i]);
      dy = Pbc(_y - y[i]);
      dz = Pbc(_z - z[i]);

      dr = sqrt(dx * dx + dy * dy + dz * dz);

      if (dr < rcut) {
        ene += 1. / pow(dr, 12) - 1. / pow(dr, 6);
      }
    }
  }
  return 4. * ene;
}

void mcmd::Move() {
  int o = 0;
  double p = 0., x_old = 0., y_old = 0., z_old = 0., x_new = 0., y_new = 0.,
         z_new = 0.;

  for (int i{}; i < npart; ++i) {
    o = static_cast<int>(rnd.Rannyu() * npart);

    x_old = x[o];
    y_old = y[o];
    z_old = z[o];

    x_new = Pbc(x[o] + delta * (rnd.Rannyu() - 0.5));
    y_new = Pbc(y[o] + delta * (rnd.Rannyu() - 0.5));
    z_new = Pbc(z[o] + delta * (rnd.Rannyu() - 0.5));
  }

  p = exp(beta + (Boltzmann(x_old, y_old, z_old, o) -
                  Boltzmann(x_new, y_new, z_new, o)));
  if (p >= rnd.Rannyu()) {
    x[o] = x_new;
    y[o] = y_new;
    z[o] = z_new;

    accepted++;
  }
  attempted++;
}
void mcmd::Measure(unsigned int istep) {
  double v = 0., w = 0., vij = 0., wij = 0., dx = 0., dy = 0., dz = 0., dr = 0.;

  for (auto &el : g_histo)
    el = 0.;

  for (int i{}; i < npart - 1; ++i) {
    for (int j{1}; j < npart; ++j) {
      // distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = sqrt(dx * dx + dy * dy + dz * dz);

      for (int k{}; k < nbins; ++k) {
        if (dr > bin_size * k && dr < bin_size * (k + 1)) {
          g_histo[k] = 3. / (2 * M_PI);
        }
      }
    }
  }
}