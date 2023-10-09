/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of bessel_i0.                                           *
 *                                                                            *
 *  bessel_i0 is free software: you can redistribute it and/or modify         *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  bessel_i0 is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with bessel_i0.  If not, see <https://www.gnu.org/licenses/>.       *
 ******************************************************************************/
#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_special_functions_real.h>

/*  Coefficients for the Remez polynomial.                                    */
#define A0 (3.98942651588301770E-01F)
#define A1 (4.98327234176892844E-02F)
#define A2 (2.91866904423115499E-02F)
#define A3 (1.35614940793742178E-02F)
#define A4 (1.31409251787866793E-01F)

#define TMPL_POLY_EVAL(z) A0 + z*(A1 + z*(A2 + z*(A3 + z*A4)))

TMPL_INLINE_DECL
float tmpl_Float_Bessel_I0_Holoborodko(float x)
{
    const float rcpr_x = 1.0F / x;
    const float exp_x = tmpl_Float_Exp_Pos_Kernel(x);
    const float inv_sqrt_x = tmpl_Float_Sqrt(rcpr_x);
    const float poly = TMPL_POLY_EVAL(rcpr_x);
    return exp_x * poly * inv_sqrt_x;
}
/*  End of tmpl_Float_Bessel_I0_Holoborodko.                                  */

#undef A0
#undef A1
#undef A2
#undef A3
#undef A4
#undef TMPL_POLY_EVAL
