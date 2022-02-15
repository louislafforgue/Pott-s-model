class Lattice {
  public:
    Lattice(int,int,double, double [], double []);      // Constructor
    int get_spin_index(int, int);
    double get_spin_value(double [], int, int);
    int get_right(int);
    int get_left(int);
    int get_above(int);
    int get_below(int);
    double get_energy(double []);
    double get_magnetization(double []);
    void Metropolis(double [], double, int);
    
  private:
    int L;
    int q;
    double h;  
};
