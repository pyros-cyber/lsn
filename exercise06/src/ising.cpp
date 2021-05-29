#include "ising.h"

Ising1D::Ising1D(string old_configuration) : rnd{"../Primes", "../seed.in"} {
  props = {"energy", "capacity", "magnetization", "susceptibility"};
  for (auto &elem : props) {
    walker[elem] = 0.;
    block_average[elem] = 0.;
    glob_average[elem] = 0.;
    glob_average2[elem] = 0.;
  }

  Input();
  /*
    if(metro) {
      Move = &Ising1D::MetropolisMove;
      std::cout << "Metropolis sampling " << std::endl;
    } else {
      Move = &Ising1D::GibbsMove;
      std::cout << "Gibbs sampling " << std::endl;
    }
  */
  if (old_configuration == "") { // start from random config
    // initial configuration
    spin_conf.resize(n_spin);
    for (auto &it : spin_conf) {
      if (rnd.Rannyu() >= 0.5)
        it = 1;
      else
        it = -1;
    }
  } else {
    ifstream input_file(old_configuration);
    if (input_file.is_open()) {
      double tmp;
      input_file >> tmp;
      while (!input_file.eof()) {
        spin_conf.push_back(tmp);
        input_file >> tmp;
      }
    } else {
      cerr << "ERROR: unable to open " << old_configuration << "\n";
      exit(1);
    }
    input_file.close();
    if (spin_conf.size() != n_spin) {
      cerr << "Spin configurations found  in '" << old_configuration
           << "' don't belong to a simulation compatible with "
           << "parameters defined in 'input.dat'\n";
      exit(0);
    }
  }

  // Evaluate energy etc. of the initial configuration
  Measure();

  // Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker["energy"] / n_spin << endl;
}

void Ising1D::Input() {
  ifstream ReadInput;
  if (ReadInput.fail()) {
    cerr << "ERROR: unable to open 'input.dat'\n";
    exit(1);
  }

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  // Read seed for random numbers
  int p1, p2;
  ifstream Primes("../Primes");
  Primes >> p1 >> p2;
  Primes.close();

  ifstream input("../seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  input.close();

  // Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0 / temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> n_spin;
  cout << "Number of spins = " << n_spin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs
  ReadInput >> nblk;
  ReadInput >> nstep;

  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();
}

void Ising1D::Measure(int step) {
  // int bin;
  double u = 0.0, m = 0.0;

  // cycle over spins
  for (unsigned int i = 0; i < n_spin; ++i) {
    u += -J * spin_conf[i] * spin_conf[Pbc(i + 1)] -
         0.5 * h * (spin_conf[i] + spin_conf[Pbc(i + 1)]);
    m += spin_conf[i];
  }
  walker.at("energy") = u;
  walker.at("capacity") = u * u;
  walker.at("magnetization") = m;
  walker.at("susceptibility") = m * m;

  if (step != 0) {
    //    std::cout << "Acceptance rate " << accepted/attempted << std::endl <<
    //    std::endl;
    ofstream EneIn, CapIn, MagIn, SusIn;
    EneIn.open("results/instant.ene.0");
    CapIn.open("results/instant.cap.0");
    MagIn.open("results/instant.mag.0");
    SusIn.open("results/instant.sus.0");

    if (EneIn.is_open() && CapIn.is_open() && MagIn.is_open() &&
        SusIn.is_open()) {
      EneIn << step << " " << walker.at("energy") << endl;
      CapIn << step << " " << walker.at("capacity") << endl;
      MagIn << step << " " << walker.at("magnetization") << endl;
      SusIn << step << " " << walker.at("susceptibility") << endl;
    } else {
      cerr << "ERROR: unable to open output files.\nExiting.\n";
      exit(1);
    }

    EneIn.close();
    CapIn.close();
    MagIn.close();
    SusIn.close();
  }
}

void Ising1D::Reset(int iblk) {
  // Reset block averages
  // we choose to lookup the maps' elements with .at for safety:
  // If you access a key using the indexing operator [] that is
  // not currently a part of a map, then it automatically adds
  // a key for you!!
  if (iblk == 0) {
    for (auto &el : props) {
      glob_average.at(el) = 0.;
      glob_average2.at(el) = 0.;
    }
  }

  for (auto &el : props)
    block_average.at(el) = 0.;

  blk_norm = 0.;
  attempted = 0.;
  accepted = 0.;
}

void Ising1D::Accumulate() {
  // Update block averages
  for (auto &el : props) {
    block_average.at(el) += walker.at(el);
  }
  blk_norm++;
}

void Ising1D::BlockAverages(int iblk) {
  // Print results for current block
  ofstream Ene, Heat, Mag, Chi;
  Ene.open("results/ene.dat"); // delete content of file
  Heat.open("results/heat.dat");
  Mag.open("results/mag.dat");
  Chi.open("results/chi.dat");

  if (Ene.is_open() && Heat.is_open() && Mag.is_open() && Chi.is_open()) {
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted / attempted << endl << endl;

    stima_u = (block_average.at("energy") / blk_norm) / n_spin;
    glob_average.at("energy") += stima_u;
    glob_average2.at("energy") += stima_u * stima_u;
    err_u = Error(glob_average.at("energy"), glob_average2.at("energy"), iblk);
    Ene << " " << iblk << " " << glob_average.at("energy") / iblk << " "
        << err_u << endl;
    Ene.close();

    stima_c = (block_average.at("capacity") / blk_norm) / n_spin;
    stima_c = (stima_c - stima_u * stima_u * n_spin) / (temp * temp);
    //  stima_c = beta*beta*((block_average.at("capacity")/blk_norm)
    //  -(std::pow(block_average.at("energy"),2)/blk_norm)); stima_c /=
    //  (double)n_spin;
    glob_average.at("capacity") += stima_c;
    glob_average2.at("capacity") += stima_c * stima_c;
    err_c =
        Error(glob_average.at("capacity"), glob_average2.at("capacity"), iblk);
    Heat << " " << iblk << " " << glob_average.at("capacity") / iblk << " "
         << err_c << endl;
    Heat.close();

    stima_m = (block_average.at("magnetization") / blk_norm) / n_spin;
    glob_average.at("magnetization") += stima_m;
    glob_average2.at("magnetization") += stima_m * stima_m;
    err_m = Error(glob_average.at("magnetization"),
                  glob_average2.at("magnetization"), iblk);
    Mag << " " << iblk << " " << glob_average.at("magnetization") / iblk << " "
        << err_m << endl;
    Mag.close();

    if (h != 0)
      stima_x = 0.;
    else
      stima_x =
          ((block_average.at("susceptibility") / blk_norm) / n_spin) / temp;

    stima_x = (beta * block_average.at("susceptibility") / blk_norm) / n_spin;
    glob_average.at("susceptibility") += stima_x;
    glob_average2.at("susceptibility") += stima_x * stima_x;
    err_x = Error(glob_average.at("susceptibility"),
                  glob_average2.at("susceptibility"), iblk);
    Chi << " " << iblk << " " << glob_average.at("susceptibility") / iblk << " "
        << err_x << endl;
    Chi.close();
  } else {
    cerr << "ERROR: unable to open output files.\nExiting.\n";
    exit(1);
  }
}

void Ising1D::ConfFinal() {
  ofstream WriteConf("config.final");

  if (WriteConf.is_open()) {
    cout << "Print final configuration to file config.final " << endl << endl;
    for (auto &el : spin_conf) {
      WriteConf << el << std::endl;
    }
  } else {
    cerr << "ERROR: unable to open output files.\nExiting.\n";
    exit(1);
  }
  WriteConf.close();
  rnd.SaveSeed();
}

void Ising1D::MetropolisMove() {
  int flip;
  double en, r;

  for (int i{}; i < n_spin; ++i) {
    // Select randomly a particle (for C++ syntax, 0 <= flip <= n_spin-1)
    flip = static_cast<int>(rnd.Rannyu() * n_spin);
    en = -2 * Boltzmann(spin_conf[flip], flip);
    en = exp(-beta * en);
    double alpha = min(1., en);
    r = rnd.Rannyu();
    if (r < alpha) {
      spin_conf[flip] = -spin_conf[flip];
      accepted++;
    }
    attempted++;
  }
}

void Ising1D::GibbsMove() {
  int flip;
  double p, en, en_up, en_down, Q, r;
  for (int i{}; i < n_spin; ++i) {
    // Select randomly a particle (for C++ syntax, 0 <= o <= n_spin-1)
    flip = static_cast<int>(rnd.Rannyu() * n_spin);
    en = Boltzmann(spin_conf[flip], flip);
    en_up = Boltzmann(1, flip);
    en_down = Boltzmann(-1, flip);
    Q = exp(-beta * en_up) + exp(-beta * en_down);
    p = exp(-beta * en) / Q;
    r = rnd.Rannyu();
    if (p > r) {
      spin_conf[flip] = 1.;
    } else {
      spin_conf[flip] = -1.;
    }
    accepted++;
    attempted++;
  }
}
/*
void Ising1D::Run(bool instant) {
  Input();                                 // Inizialization
  for (int iblk = 1; iblk <= nblk; ++iblk) // Simulation
  {
    Reset(iblk); // Reset block averages
    for (int istep = 1; istep <= nstep; ++istep) {
      (*this.*Move)();
      if (instant)
        Measure(istep);
      else
        Measure();
      Accumulate(); // Update block averages
    }
    BlockAverages(iblk); // Print results for current block
  }
  ConfFinal(); // Write final configuration
}*/