/**
 * \file MagneticCircle.cpp
 * \brief Implementation for GeographicLib::MagneticCircle class
 *
 * Copyright (c) Charles Karney (2011-2015) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#include <GeographicLib/MagneticCircle.hpp>
#include <fstream>
#include <sstream>
#include <GeographicLib/Geocentric.hpp>

namespace GeographicLib {

  using namespace std;

  void MagneticCircle::Field(real lon, bool diffp,
                             real& Bx, real& By, real& Bz,
                             real& Bxt, real& Byt, real& Bzt) const {
    real clam, slam;
    CircularEngine::cossin(lon, clam, slam);
    real M[Geocentric::dim2_];
    Geocentric::Rotation(_sphi, _cphi, slam, clam, M);
    real BX0, BY0, BZ0, BX1, BY1, BZ1; // Components in geocentric basis
    real BXc = 0, BYc = 0, BZc = 0;
    _circ0(clam, slam, BX0, BY0, BZ0);
    _circ1(clam, slam, BX1, BY1, BZ1);
    if (_constterm)
      _circ2(clam, slam, BXc, BYc, BZc);
    if (_interpolate) {
      BX1 = (BX1 - BX0) / _dt0;
      BY1 = (BY1 - BY0) / _dt0;
      BZ1 = (BZ1 - BZ0) / _dt0;
    }
    BX0 += _t1 * BX1 + BXc;
    BY0 += _t1 * BY1 + BYc;
    BZ0 += _t1 * BZ1 + BZc;
    if (diffp) {
      Geocentric::Unrotate(M, BX1, BY1, BZ1, Bxt, Byt, Bzt);
      Bxt *= - _a;
      Byt *= - _a;
      Bzt *= - _a;
    }
    Geocentric::Unrotate(M, BX0, BY0, BZ0, Bx, By, Bz);
    Bx *= - _a;
    By *= - _a;
    Bz *= - _a;
  }

} // namespace GeographicLib
