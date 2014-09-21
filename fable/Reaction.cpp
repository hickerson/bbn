/**
 * Reaction.cpp
 * 
 * Author: Kevin Peter Hickerson
 * created: Sep 12, 2014
 */
#include "Reaction.hpp"

/*
    std::cout << "reaction type: "<< type << "\n";
    int ri = s[0][type-1]; 		/// # of incoming nuclide i.
    int rj = s[1][type-1]; 		/// # of incoming nuclide j.
    int rk = s[2][type-1]; 		/// # of outgoing nuclide k.
    int rl = s[3][type-1]; 		/// # of outgoing nuclide l.
    //..........COMPUTE DIFFERENT REACTION RATES.
    switch (type) {
        case 1: /// 1-0-0-1 configuration.
            ci = f[n]; 			/// (Ref 1).
            cj = 0;
            ck = 0;
            cl = r[n];
            break;

        case 2: /// 1-1-0-1 configuration.
            r[n] = rev(n) * 1.e+10f * T932 * ex(-q9(n) / T9) * f[n]; 	/// (Ref 2).
            f[n] = rhob * f[n];
            ci = y(j) * f[n] / 2.;
            cj = y(i) * f[n] / 2.;
            ck = 0;
            cl = r[n];
            break;

        case 3: /// 1-1-1-1 configuration.
            f[n] = rhob * f[n]; 			/// (Ref 3).
            r[n] = rev(n) * ex(-q9(n) / T9) * f[n];
            ci = y(j) * f[n] / 2;
            cj = y(i) * f[n] / 2;
            ck = y(l) * r[n] / 2;
            cl = y(k) * r[n] / 2;
            break;

        case 4: /// 1-0-0-2 configuration.
            ci = f[n];
            cj = 0;
            ck = 0;
            cl = y(l) * r[n] / 2;
            break;

        case 5: /// 1-1-0-2 configuration.
            f[n] = rhob * f[n];
            r[n] = rev(n) * ex(-q9(n) / T9) * f[n]; 	/// (Ref 3).
            ci = y(j) * f[n] / 2;
            cj = y(i) * f[n] / 2;
            ck = 0;
            cl = y(l) * r[n] / 2;
            break;

        case 6: /// 2-0-1-1 configuration.
            f[n] = rhob * f[n];
            r[n] = rev(n) * ex(-q9(n) / T9) * f[n]; 	/// (Ref 3).
            ci = y(i) * f[n] / 2;
            cj = 0;
            ck = y(l) * r[n] / 2;
            cl = y(k) * r[n] / 2;
            break;

        case 7: /// 3-0-0-1 configuration.
            //(Ref 4).
            r[n] = rev(n) * 1.e+20f * T932 * T932 * ex(-q9(n) / T9) * f[n];
            f[n] = rhob * rhob * f[n];
            ci = y(i) * y(i) * f[n] / 6;
            cj = 0;
            ck = 0;
            cl = r[n];
            break;

        case 8: /// 2-1-0-1 configuration.
            //(Ref 4).
            r[n] = rev(n) * 1.e+20f * T932 * T932 * ex(-q9(n) / T9) * f[n];
            f[n] = rhob * rhob * f[n];
            ci = y(j) * y(i) * f[n] / 3.;
            cj = y(i) * y(i) * f[n] / 6.;
            ck = 0.;
            cl = r[n];
            break;

        case 9: /// 1-1-1-2 configuration.
            f[n] = rhob * f[n];
            //(Ref 5)
            r[n] = rev(n) * 1.e-10f * T9m32 * rhob * ex(-q9(n) / T9) * f[n];
            ci = y(j) * f[n] / 2.;
            cj = y(i) * f[n] / 2.;
            ck = y(l) * y(l) * r[n] / 6.;
            cl = y(k) * y(l) * r[n] / 3.;
            break;

        case 10: /// 1-1-0-3 configuration.
            f[n] = rhob * f[n];
            //(Ref 5)
            r[n] = rev(n) * 1.e-10f * T9m32 * rhob * ex(-q9(n) / T9) * f[n];
            ci = y(j) * f[n] / 2.;
            cj = y(i) * f[n] / 2.;
            ck = 0.;
            cl = y(l) * y(l) * r[n] / 6.;
            break;

        case 11: /// 2-0-2-1 configuration.
            f[n] = rhob * f[n];
            //(Ref 5)
            r[n] = rev(n) * 1.e-10f * T9m32 * rhob * ex(-q9(n) / T9) * f[n];
            ci = y(i) * f[n] / 2.;
            cj = 0.;
            ck = y(l) * y(k) * r[n] / 3.;
            cl = y(k) * y(k) * r[n] / 6.;
            break;

        default: break;
    }
*/
