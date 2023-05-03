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

#ifndef VISPROPT_H
#define VISPROPT_H

#include <QWidget>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>

using namespace cv;

namespace Ui {
class VISproPT;
}

class VISproPT : public QWidget
{
    Q_OBJECT

public:
    explicit VISproPT(QWidget *parent = nullptr);
    ~VISproPT();
    void setWin(const std::string& _winname);

private:
    Ui::VISproPT *ui;
    void on_mouse_internal(int ev, int x, int y);
    std::string winname;
    friend void on_mouse(int ev, int x, int y, int, void* obj);

public slots:
    void setCamera();
    void selectImg();
    void viewImg(QString fileImgCam1,QString fileIngCam2);
    void setIntensify();
    void process();
    void grabFrame();
    void GoTo();
    void scan();
    void saveAcam(QString file2save);
    void readAcam(QString file2read);
    void callReadAcam();
    void callSaveAcam();
    void stat(int iq);
    void intFat();
    void expo(QString file2expo,int iq);
    void doIntFat();
    void CALLcalib1();
    void CALLcalib2();
    void calibrate(int iCam);
    void res2GUI();
    void idealShape();
    void idealHight(int iWarning);
    void shape();
    void slope();
    void rmTilt();
    void mapMat();
    double slopeChi2();
    void readParam();
    void setFrmCam1();
    void setFrmCam2();
    void inpainting();
    void shapeComputing();
    void compSou();
    void setCB();
    void refBFM();
};

#endif // VISPROPT_H
