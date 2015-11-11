#ifndef NUMDEV_H
#define NUMDEV_H


class numdev{
    public:
        numdev();
        virtual ~numdev();
        static double lap5p_o2(double* f, int m, double h);
        static double jac9p_o2(double* p, double* z, int m, double h);


};

#endif
