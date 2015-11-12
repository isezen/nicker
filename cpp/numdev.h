#ifndef NUMDEV_H
#define NUMDEV_H


class numdev{
    public:
        numdev();
        virtual ~numdev();
        static void lap5p_o2(double* a, double* t, int l, int m, double h);
        static void jac9p_o2(double* a, double* p, double* z, int l, int m, double h);
        static void RELAX1(double* p, double* eta, int n, int m, double alpha);

    private:
        static void j1_first_col (double* a, double* p, double* z, int l, int m, double h);
        static void j1_last_col (double* a, double* p, double* z, int l, int m, double h);
        static void j1 (double* a, double* p, double* z, int l, int m, double h);
};

#endif
