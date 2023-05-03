/*Author: Marco Montecchi
          ENEA-Casaccia
          Roma, Italy
          email: marco.montecchi@enea.it

This software covers all the steps of the 3D shape measurement of parabolic-trough reflective panels acchomplished with the innovative instrument named VISproPT:
    1) Camera-lens calibration, for image undistortion
    2) Instrument calibration
    3) Image-processing for evaluating: i) 3D shape (slopes dz/dx, dz/dy and height z), ii) deviations from the ideal shape and, last but not least, iii) evaluation of the intercept factor at a given longitudinal angle.

   Copyright (C) 2023  Marco Montecchi

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include "vispropt.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    VISproPT w;
    w.show();

    return a.exec();
}
