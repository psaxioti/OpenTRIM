#ifndef _ELEMENTS_H_
#define _ELEMENTS_H_

/* Extracted from NIST data in May 2007 (physics.nist.gov/PhysRefData/Compositions/index.html)
Developers and Contributors:
J. S. Coursey, D. J. Schwab, and R. A. Dragoset
NIST, Physics Laboratory, Office of Electronic Commerce in Scientific and Engineering Data
(There are 100 data but only 92 are used with SRIM 2006) */

struct elements {
    static const int max_atomic_num = 92;
    static float mostAbundantIsotope(int atomic_num);
    static const char* name(int atomic_num);
    static float mass(int atomic_num);
};

#endif // ELEMENTS_H
