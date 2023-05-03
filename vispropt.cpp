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
#include "ui_vispropt.h"
#include <QtGui>
#include <QFile>
#include <QFileDialog>
#include <QMessageBox>
#include <QInputDialog>
#include <QTextStream>
#include <QTime>
#include <QDate>
#include <fstream>
#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <cminpack.h>
#include <random>

using namespace std;
using namespace cv;


//global variables *******************************
QString dir="/run/media/marco/d2quadra/VISproPT";//su Z440
//QString dir="/home/marco/Workspace/VISproPT"; //su portatile
//QString dir="/home/operatore/Documenti/VISproPT"; //su HP_xw4600
//QString dir="/home/hpz210/Workspace/VISproPT";//HP Z210

QString filetXYZtargets=dir+"/xyzTargets.txt";// <<-targets2 =black/white
QString fileAcam=dir+"/AcamNew.txt";
QString fileMatDz=dir+"/expo/MatrixDz.txt";
QString fileMatRMSz=dir+"/expo/MatrixZrms.txt";
QString fileMatXslope=dir+"/expo/MatrixDslopeX.txt";
QString fileMatYslope=dir+"/expo/MatrixDslopeY.txt";
QString fileMatIntFat=dir+"/expo/MatrixIntFat.txt";
QString PathSave=dir+"/frames/";
QString PathPro="";
QString part1cmd="ssh root@192.168.125.130 bash /retinae/prj/cmd/script/";
QString filename1,filename2;
QString correz=dir+"/Correzioni.txt";
QString dirMeas;
int mpxX1,mpxY1,mpxX2,mpxY2,drag;
int iReading=0;
int iXcur,iXlab;
int iWr=0;//0->no verbose
double XlabCam1,XlabCam2,x4set;
double px2mm;//pixel -> mm
double pig=acos(-1);
double chi2fin,chi2min;
double Vservice[3];//vettore di servizio
double Mean,standDev;//for statistics
double MinMax[3][2];
double minDis;

//correzione artefatti
int iArtCor;
double corArte[216][5];
/*
    corArte[i][0]=xMotor (mm)
    corArte[i][1]=dZcam1 (mm)
    corArte[i][2]=dZcam2 (mm)
    corArte[i][3]=dPitch (rad)
    corArte[i][4]=dRoll  (rad)
*/

//result matrix where store the data of each one of the observed reflected points
double results[550000];
/*  results[i*6+ 0] = N. frame
    results[i*6+ 1] = N. Camera
    results[i*6+ 2] = N. source
    results[i*6+ 3] =
    results[i*6+ 4] = j bubble-coordinate (px)
    results[i*6+ 5] = i bubble-coordinate (px) */
int iTot=0;//N. data in result matrix

//sampling matrix
const int Ni=301;
const int Nj=201;
int nScellFill=0;
double S[Nj][Ni][12]={};
/*  S[j][i][0] = vn[0] \
    S[j][i][1] = vn[1] - normal unit vector to the surface
    S[j][i][2] = vn[2] /
    S[j][i][3] = z_ideal (height) (mm)
    S[j][i][4] = N. of independent evaluations
    S[j][i][5] = dz/dx (partial y-derivative) // i
    S[j][i][6] = dz/dy (partial x-derivative) // j
    S[j][i][7] = z_exp (height) (mm)
    S[j][i][8] = dz/dx_ideal
    S[j][i][9] = dz/dy_ideal
    S[j][i][10]= intFat
    S[j][i][11]= z_exp rms   */
double DS=10.0;//(mm) size of S cell
const int iScenter=(Ni-1)/2;
const int jScenter=(Nj-1)/2;
int isMIN=0;
int isMAX=Ni-1;
int jsMIN=0;
int jsMAX=Nj-1;

//attaching points & panel parameters
const int jAP=2;// N. attaching points along j axis
const int iAP=2;// N. attaching points along i axis
int jiAP[jAP][iAP][3]={};//(j,i) intex of the attaching points
//jiAP[jAP][iAP][0] = j attaching point (jap,iap)
//jiAP[jAP][iAP][1] = i attaching point (jap,iap)
//jiAP[jAP][iAP][2] = 1 if attaching point (jap,iap) is in the sampled region or active, else is 0
double Djiz[3];//largest hight-difference on the attaching points
double focal=1810.;// focal lenght
double alpha;//angle between LAB and parabola reference-frames

// BAUMER  parameters **************************************************************************
int Width= 2448;
int Height=2048;
double expTime=15000;
double fcam_1=8.0; //camera focal (mm)
double fcam_2=8.0; //camera focal (mm)
int NpxY=2048;//N. px lato Y                   (0,0) ---> x' (j)
int NpxX=2448;//N. px lato X                     |
double ccdY=7.0656;//lunghezza lato Y ccd        |
double ccdX=8.4456;//lunghezza lato X ccd        V y' (i)
double pxdimX=ccdX/double(NpxX);
double pxdimY=ccdY/double(NpxY);
double fx_1=fcam_1/pxdimX;
double fy_1=fcam_1/pxdimY;
double cx_1=NpxX/2.;
double cy_1=NpxY/2.;
double fx_2=fcam_2/pxdimX;
double fy_2=fcam_2/pxdimY;
double cx_2=NpxX/2.;
double cy_2=NpxY/2.;

Mat cameraMatrix_1 =  (Mat_<double>(3,3) <<  fx_1,          0 ,     cx_1,
                            0,fy_1 ,     cy_1,
                            0,          0 ,      1);

Mat distCoeffs_1 = (Mat_<double>(14,1)   <<    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);

Mat cameraMatrix_2 =  (Mat_<double>(3,3) <<  fx_2,          0 ,     cx_2,
                            0,fy_2 ,     cy_2,
                            0,          0 ,      1);

Mat distCoeffs_2 = (Mat_<double>(14,1)   <<    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
//***************************************************************************************************

int setFrameCam[2][2][2];
//             [0][0][iCam-1]=j_left_up_corner
//             [0][1][iCam-1]=i_left_up_corner
//             [1][0][iCam-1]=DELTAi
//             [1][1][iCam-1]=DELTAj
Point2i Pt1,Pt2;
Mat img,tmpImg;// image
Mat img8b(Height,Width,CV_8UC1);
Mat img8RGB(Height,Width,CV_8UC3);
//rgb2gray 0.299*R+0.587*G+0.114*B
//Note: OpenCV uses BGR color order, i.e. [0]=B [1]=G [2]=R
int imgN=0;
int Lheight,Lwidth;
int undst1,undst2;
double Dz[2][2];
cv::SimpleBlobDetector::Params params;

double pxX,pxY;//coordinate pixel j,i
// cameras
double Pc_1[3];//coordinate camera1
double EA_1[3];//yaw pitch roll camera1
double Pc_2[3];//coordinte monitor2
double EA_2[3];//yaw pitch roll camera2
double Dx21,Dz21;
//double Xoffset_Pc_2=-3.2;//da misure Leica 7 Giugno 2022
//double Zoffset_Pc_2=1.6; //idem
//WARNING: yaw, pitch, roll must be evaluated to rotate the Plane to return to Earth !!!!!
int Ncorn;
int NbCam1;
int NbCam2;
double Bcam1[200][3];
double Bcam2[200][3];

const int NjZlcd=16;
const int NiZlcd=11;
double MatZlcd[NjZlcd][NiZlcd][3];

//point source array
double Ps[3];//coordinate centre source
double Es[3];//yaw pitch roll source
double source[200][3];//x,y,z coordinates of point sources
double zMirror=0;//of mirror used to reflect point source in setCam

//targetsFptf
double P[2000][3];//coordinate x,y,z targets
double px[3000][4]; //coordinate in pixel dei target nell'immagine e del punto calcolato

//scanPar
double Xmin,Xmax,Xstep,Xstart;

//solarVector
double uvS[3];

//focal point
double pF[3];

//invoked functions ********************************
int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag);
int fcnRefCam(void *p, int m, int n, const double *x, double *fvec, int iflag);
void kernelpro(int icam);
int kernelSetCam(int iSetCamMeth);
void InfoFit(int info);
//void pxpy2xyz(double di,double dj);
void Plane2Earth(double yaw,double pitch, double roll,double xd, double yd, double zd);
void Earth2Plane(double yaw,double pitch, double roll,double xd, double yd, double zd);
void unitV12(double x1, double y1, double z1, double x2, double y2, double z2);
void sumV12(double x1, double y1, double z1, double x2, double y2, double z2);
void xyz2pxpy(int iCam, double X, double Y, double Z);
void PxPyDistorted(int iCam,double xTrue, double yTrue);
void pxColor(double val);
void on_mouse(int ev, int x, int y, int, void* obj){
    VISproPT* app = static_cast<VISproPT*>(obj);
    if (app)
        app->on_mouse_internal(ev, x, y);
}
void slopeComputing();
int getPosition();
int NINT(double x);

struct pointToFit2{
    int Nt;
}pTF2[1];



VISproPT::VISproPT(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::VISproPT)
{
    ui->setupUi(this);
    //signals & slots
    connect( ui->pB_Load,        SIGNAL( clicked() ),           this,SLOT(selectImg()   ));
    connect( ui->pB_aquire,      SIGNAL( clicked() ),           this,SLOT(grabFrame()   ));
    connect( ui->sB_gL,          SIGNAL(valueChanged(int)),     this,SLOT(setIntensify()));
    connect( ui->checkBox_undist,SIGNAL(stateChanged(int)),   this,SLOT(setIntensify()));
    connect( ui->pB_setCamSou,   SIGNAL( clicked() ),           this,SLOT(setCamera()   ));
    connect( ui->pB_goto,        SIGNAL( clicked() ),           this,SLOT(GoTo()        ));
    connect( ui->pB_pro,         SIGNAL( clicked() ),           this,SLOT(process()     ));
    connect( ui->pB_replotCMap,  SIGNAL( clicked() ),           this,SLOT(mapMat()      ));
    connect( ui->pB_scan,        SIGNAL(clicked(bool)),         this,SLOT(scan()        ));
    connect( ui->pB_saveAcam,    SIGNAL(clicked(bool)),         this,SLOT(callSaveAcam()));
    connect( ui->pB_calibrate_1, SIGNAL(clicked(bool)),         this,SLOT(CALLcalib1()  ));
    connect( ui->pB_calibrate_2, SIGNAL(clicked(bool)),         this,SLOT(CALLcalib2()  ));
    connect( ui->pB_idealShape,  SIGNAL(clicked(bool)),         this,SLOT(idealShape()  ));
    connect( ui->pB_shape,       SIGNAL(clicked(bool)),         this,SLOT(shape()       ));
    connect( ui->pB_slope,       SIGNAL(clicked(bool)),         this,SLOT(slope()       ));
    connect( ui->pB_rmTilt,      SIGNAL(clicked(bool)),         this,SLOT(rmTilt()      ));
    connect(ui->pB_readAcam,     SIGNAL(clicked(bool)),         this,SLOT(callReadAcam()));
    connect(ui->pushButton_setFrmCam1,SIGNAL(clicked()),        this,SLOT(setFrmCam1()  ));
    connect(ui->pushButton_setFrmCam2,SIGNAL(clicked()),        this,SLOT(setFrmCam2()  ));
    connect(ui->pushButton_compSou,SIGNAL(clicked()),           this,SLOT(compSou()     ));
    connect(ui->comboBox_chessType,SIGNAL(currentIndexChanged(int)),this,SLOT(setCB()   ));
    connect(ui->pushButton_refBFM,SIGNAL(clicked(bool)),       this,SLOT(refBFM()));

    printf("*******************************************\n");
    printf("VISproPT software\nAuthor: Marco Montecchi\n\tENEA-Casaccia, Roma, Italy\n\temail: marco.montecchi@enea.it\nVersion 3 May 2023\n");
    printf("License: The code source files are distributed as open source software \n\t under the GNU General Public License \n\t as published by the Free Software Foundation version 3\n");
    printf("*******************************************\n");

    //basic initialization
    waitKey(100);
    callReadAcam();
    ui->progressBar->setValue(0);

    QFile coAr(correz);
    if(!coAr.open(QIODevice::ReadOnly | QIODevice::Text))
        printf("file %s non trovato!!!\n",correz.toStdString().c_str());
    else{
        QTextStream stream(&coAr);
        QString line=stream.readLine();
        for(int i=0;i<216;i++){
            stream>>corArte[i][0]>>corArte[i][1]>>corArte[i][2]>>corArte[i][3]>>corArte[i][4];
            //printf("%f\t%f\t%f\t%f\t%f\n",
            //       corArte[i][0],corArte[i][1],corArte[i][2],corArte[i][3],corArte[i][4]);
        }
    }

    //get x position
    iXcur=getPosition();
    printf("iXcur= %d\n",iXcur);
    ui->lineEdit_Xcur->setText(QString::number(iXcur));

    //windows
    namedWindow("Camera1",WINDOW_NORMAL );
    namedWindow("Camera2",WINDOW_NORMAL );
    namedWindow("Board",WINDOW_NORMAL);
    namedWindow("ROIcam1",WINDOW_NORMAL);
    namedWindow("ROIcam2",WINDOW_NORMAL);
    filename1=dir+"/imgCam1.bmp";
    filename2=dir+"/imgCam2.bmp";
    Mat img1=imread((dir+"/imgCam1.bmp").toStdString(),IMREAD_ANYDEPTH);
    Mat img2=imread((dir+"/imgCam2.bmp").toStdString(),IMREAD_ANYDEPTH);
    Height=img1.rows;
    Width=img1.cols;
    printf("img: height=%d width=%d\n",Height,Width);
    imshow("Camera1",img1);
    imshow("Camera2",img2);
    setWin("Camera1");//abilita il cursore sulla finestra
    setWin("Camera2");
    setWin("ROIcam1");
    setWin("ROIcam2");
    setCB();
}

VISproPT::~VISproPT()
{
    saveAcam(fileAcam);
    destroyAllWindows();
    delete ui;
}

void VISproPT::compSou(){
    //setting of point source coordinates
    double Dsource=ui->dSB_pSourceStep->value();//point source step
    for(int i=0;i<200;i++){
        for(int k=0;k<3;k++)
            source[i][k]=0.;
    }
    source[99][0]=-Dsource;                       //         ^ x
    source[100][0]=Dsource;                       //         |
    for(int i=0;i<=98;i++){                       //         |
        source[98-i][1]=-Dsource*(0.5+double(i)); //         |
        source[101+i][1]=Dsource*(0.5+double(i)); //  <------0
    }                                             //  y
    Ps[0]= ui-> dSB_Xs -> value();
    Ps[1]= ui-> dSB_Ys -> value();
    Ps[2]= ui-> dSB_Zs -> value();
    Es[0]= ui-> dSB_Yaw_s -> value(); //yaw
    Es[1]= ui-> dSB_Pitch_s -> value(); //pitch
    Es[2]= ui-> dSB_Roll_s -> value(); //roll
    for(int i=0;i<3;i++)
        Es[i]= Es[i]/180.*pig;
    int nSource;
    for(int k=0;k<4;k++){
        if(k==0) nSource=5;
        else if(k==1) nSource=99;
        else if(k==2) nSource=100;
        else if(k==3) nSource=193;
        Plane2Earth(Es[0],Es[1],Es[2],
            source[nSource][0],source[nSource][1],source[nSource][2]);
        printf("XYZ Source[%d][k]= %f\t%f\t%f\n",
           nSource,Vservice[0]+Ps[0],Vservice[1]+Ps[1],Vservice[2]+Ps[2]);
    }
}

void VISproPT::callSaveAcam(){
    QString file2save;
    QString AcamName=ui->lineEdit_AcamName->text();
    if(AcamName.isEmpty()){
        file2save=QFileDialog::getSaveFileName(
                this,
                "Filename to save",
                dir,
                "file (*.txt)");
        AcamName=file2save.section('/',-1,-1);
        ui->lineEdit_AcamName->setText(AcamName);
    }
    else
        file2save=dir+"/"+AcamName;
    saveAcam(file2save);
}


void VISproPT::saveAcam(QString file2save){
    if(iReading==1)
        return;
    readParam();
    printf("saveAcam as %s\n",file2save.toStdString().c_str());
    Pc_1[0]= ui-> dSB_Xc_1 -> value();
    Pc_1[1]= ui-> dSB_Yc_1 -> value();
    Pc_1[2]= ui-> dSB_Zc_1 -> value();
    EA_1[0]= ui-> dSB_Yaw_1 -> value(); //yaw
    EA_1[1]= ui-> dSB_Pitch_1 -> value(); //pitch
    EA_1[2]= ui-> dSB_Roll_1 -> value(); //roll
    fcam_1 = ui-> dSB_fcam_1 ->value();//fcamera
    Pc_2[0]= ui-> dSB_Xc_2 -> value();
    Pc_2[1]= ui-> dSB_Yc_2 -> value();
    Pc_2[2]= ui-> dSB_Zc_2 -> value();
    EA_2[0]= ui-> dSB_Yaw_2 -> value(); //yaw
    EA_2[1]= ui-> dSB_Pitch_2 -> value(); //pitch
    EA_2[2]= ui-> dSB_Roll_2 -> value(); //roll
    fcam_2 = ui-> dSB_fcam_2 ->value();//fcamera
    for(int i=0;i<3;i++){
        EA_1[i]= EA_1[i]/180.*pig;
        EA_2[i]= EA_2[i]/180.*pig;
    }
    int thres= ui->sB_gL -> value();   //threshold
    expTime = ui->dSBexpTime -> value();
    Ps[0]= ui-> dSB_Xs -> value();
    Ps[1]= ui-> dSB_Ys -> value();
    Ps[2]= ui-> dSB_Zs -> value();
    Es[0]= ui-> dSB_Yaw_s -> value(); //yaw
    Es[1]= ui-> dSB_Pitch_s -> value(); //pitch
    Es[2]= ui-> dSB_Roll_s -> value(); //roll
    for(int i=0;i<3;i++)
        Es[i]= Es[i]/180.*pig;
    focal= ui-> dSB_focal -> value();
    int iPanel= ui -> comboBox_panel ->currentIndex();
    double Zvalue;
    if(iPanel==4)
        Zvalue=ui->dSB_zMirror -> value();
    else
        Zvalue=ui->doubleSpinBox_dz0hold-> value();
    Xmin=ui->dSB_Xmin->value();
    Xmax=ui->dSB_Xmax->value();
    Xstep=ui->dSB_Xstep->value();

    QFile file(file2save);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream streams ( &file );
    streams<<"Pc&EA&fcam_1\t"<<Pc_1[0]<<"\t"<<Pc_1[1]<<"\t"<<Pc_1[2]<<"\t"
          <<EA_1[0]*180./pig<<"\t"<<EA_1[1]*180./pig<<"\t"<<EA_1[2]*180./pig<<"\t"<<fcam_1<<"\n";
    streams<<"Pc&EA&fcam_2\t"<<Pc_2[0]<<"\t"<<Pc_2[1]<<"\t"<<Pc_2[2]<<"\t"
          <<EA_2[0]*180./pig<<"\t"<<EA_2[1]*180./pig<<"\t"<<EA_2[2]*180./pig<<"\t"<<fcam_2<<"\n";
    for(int i=0;i<3;i++){
        streams<<"camMat_1#\t"<<cameraMatrix_1.at<double>(i,0)<<"\t"<<cameraMatrix_1.at<double>(i,1)
              <<"\t"<<cameraMatrix_1.at<double>(i,2)<<"\n";
    }
    streams<<"distCoeff_1";
    for(int j=0;j<14;j++)
        streams<<"\t"<<distCoeffs_1.at<double>(j);
    streams<<"\n";
    for(int i=0;i<3;i++){
        streams<<"camMat_2#\t"<<cameraMatrix_2.at<double>(i,0)<<"\t"<<cameraMatrix_2.at<double>(i,1)
              <<"\t"<<cameraMatrix_2.at<double>(i,2)<<"\n";
    }
    streams<<"distCoeff_2";
    for(int j=0;j<14;j++)
        streams<<"\t"<<distCoeffs_2.at<double>(j);
    streams<<"\n";
    streams<<"greyThres\t"<<thres<<"\n";
    streams<<"ExpTime\t"<<expTime<<"\n";
    streams<<"Ps&EsPS\t"<<Ps[0]<<"\t"<<Ps[1]<<"\t"<<Ps[2]<<"\t"<<Es[0]*180./pig<<"\t"<<Es[1]*180./pig<<"\t"<<Es[2]*180./pig<<"\n";
    streams<<"panelFocal\t"<<focal<<"\n";
    streams<<"panelType\t"<<iPanel<<"\n";
    streams<<"zValue\t"<<Zvalue<<"\n";
    streams<<"xMinMaxStep\t"<<Xmin<<"\t"<<Xmax<<"\t"<<Xstep<<"\n";
    x4set=ui->dSB_x4set->value();
    streams<<"x4set\t"<<x4set<<"\n";

    file.close();
}

void VISproPT::callReadAcam(){
    QString file2read;
    QString AcamName=ui->lineEdit_AcamName->text();
    if(AcamName.isEmpty())
        file2read="";
    else
        file2read=dir+"/"+AcamName;
    readAcam(file2read);
}


void VISproPT::readAcam(QString file2read){
    QMessageBox msgBox1;
    int ret;
    if(!file2read.isEmpty()){
        msgBox1.setText("ATTENTION!");
        msgBox1.setInformativeText("Do you want to read Acam @ "+file2read+" ?");
        msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox1.setDefaultButton(QMessageBox::Yes);
        ret = msgBox1.exec();
    }
    //QString dir2;
    if(file2read.isEmpty() || ret==QMessageBox::No){
        //dir2="";
        //QString dir2=QFileDialog::getExistingDirectory(this, tr("Open Directory where read Acam"),
        //                                       dir,
        //                                       QFileDialog::ShowDirsOnly);
        QString dir2=QFileDialog::getOpenFileName(
                    this,
                    "Choose a Acam file", //titolo della finestra
                    dir, //directory iniziale
                    "file (*.txt)"); //tipi di file da cercare
        if(dir2.isEmpty()){
            printf("abort!\n");
            return;
        }
        file2read=dir2;//+"/AcamNew.txt";
    }
    iReading=1;
    QString comment;
    QFile file(file2read);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    fileAcam=file2read;
    printf("readAcam @ %s\n",file2read.toStdString().c_str());
    QString AcamName=file2read.section('/',-1,-1);
    ui->lineEdit_AcamName->setText(AcamName);
    QTextStream streamr ( &file );
    streamr>>comment>>Pc_1[0]>>Pc_1[1]>>Pc_1[2]>>EA_1[0]>>EA_1[1]>>EA_1[2]>>fcam_1;
    streamr>>comment>>Pc_2[0]>>Pc_2[1]>>Pc_2[2]>>EA_2[0]>>EA_2[1]>>EA_2[2]>>fcam_2;
    //camera1 attitude & focal length
    ui-> dSB_Xc_1    -> setValue(Pc_1[0]);
    ui-> dSB_Yc_1    -> setValue(Pc_1[1]);
    ui-> dSB_Zc_1    -> setValue(Pc_1[2]);
    ui-> dSB_Yaw_1   -> setValue(EA_1[0]); //yaw
    ui-> dSB_Pitch_1 -> setValue(EA_1[1]); //pitch
    ui-> dSB_Roll_1  -> setValue(EA_1[2]); //roll
    ui-> dSB_fcam_1  -> setValue(fcam_1);
    //camera2 attitude & focal length
    ui-> dSB_Xc_2    -> setValue(Pc_2[0]);
    ui-> dSB_Yc_2    -> setValue(Pc_2[1]);
    ui-> dSB_Zc_2    -> setValue(Pc_2[2]);
    ui-> dSB_Yaw_2   -> setValue(EA_2[0]); //yaw
    ui-> dSB_Pitch_2 -> setValue(EA_2[1]); //pitch
    ui-> dSB_Roll_2  -> setValue(EA_2[2]); //roll
    ui-> dSB_fcam_2  -> setValue(fcam_2);
    ui-> dSB_Dx21    -> setValue(Pc_2[0]-Pc_1[0]);
    ui-> dSB_Dz21    -> setValue(Pc_2[2]-Pc_1[2]);
    for(int i=0;i<3;i++){
        EA_1[i]= EA_1[i]/180.*pig;
        EA_2[i]= EA_2[i]/180.*pig;
    }
    for(int i=0;i<3;i++)
        streamr>>comment>>cameraMatrix_1.at<double>(i,0)>>cameraMatrix_1.at<double>(i,1)
                                                     >>cameraMatrix_1.at<double>(i,2);
    QString linea, st1;
    linea=streamr.readLine();
    linea=streamr.readLine();//boh!!!!!
    cout<<linea.toStdString().c_str()<<"\n";
    QStringList List;
    List=linea.split("\t");
    int nV=List.count();
    printf("QList cam1: nV=%d\n",nV);
    for(int iv=1;iv<=14;iv++){
        if(iv<nV){
            st1=List.at(iv).toLocal8Bit().constData();
            distCoeffs_1.at<double>(iv-1)=st1.toDouble();
        }
        else
            distCoeffs_1.at<double>(iv-1)=0.;
    }
    ui->lineEdit_fx_1->setText(QString::number(cameraMatrix_1.at<double>(0,0)));
    ui->lineEdit_fy_1->setText(QString::number(cameraMatrix_1.at<double>(1,1)));
    ui->lineEdit_cx_1->setText(QString::number(cameraMatrix_1.at<double>(0,2)));
    ui->lineEdit_cy_1->setText(QString::number(cameraMatrix_1.at<double>(1,2)));
    ui->lineEdit_k1_1->setText(QString::number(distCoeffs_1.at<double>(0)));
    ui->lineEdit_k2_1->setText(QString::number(distCoeffs_1.at<double>(1)));
    ui->lineEdit_p1_1->setText(QString::number(distCoeffs_1.at<double>(2)));
    ui->lineEdit_p2_1->setText(QString::number(distCoeffs_1.at<double>(3)));
    ui->lineEdit_k3_1->setText(QString::number(distCoeffs_1.at<double>(4)));
    ui->lineEdit_k4_1->setText(QString::number(distCoeffs_1.at<double>(5)));
    ui->lineEdit_k5_1->setText(QString::number(distCoeffs_1.at<double>(6)));
    ui->lineEdit_k6_1->setText(QString::number(distCoeffs_1.at<double>(7)));
    ui->lineEdit_s1_1->setText(QString::number(distCoeffs_1.at<double>(8)));
    ui->lineEdit_s2_1->setText(QString::number(distCoeffs_1.at<double>(9)));
    ui->lineEdit_s3_1->setText(QString::number(distCoeffs_1.at<double>(10)));
    ui->lineEdit_s4_1->setText(QString::number(distCoeffs_1.at<double>(11)));
    ui->lineEdit_tx_1->setText(QString::number(distCoeffs_1.at<double>(12)));
    ui->lineEdit_ty_1->setText(QString::number(distCoeffs_1.at<double>(13)));
    fx_1=cameraMatrix_1.at<double>(0,0);
    fy_1=cameraMatrix_1.at<double>(1,1);
    cx_1=cameraMatrix_1.at<double>(0,2);
    cy_1=cameraMatrix_1.at<double>(1,2);
    for(int i=0;i<3;i++)
        streamr>>comment>>cameraMatrix_2.at<double>(i,0)>>cameraMatrix_2.at<double>(i,1)
                                                     >>cameraMatrix_2.at<double>(i,2);
    linea=streamr.readLine();
    linea=streamr.readLine();
    cout<<linea.toStdString().c_str()<<"\n";
    List=linea.split("\t");
    nV=List.count();
    printf("QList cam2: nV=%d\n",nV);
    for(int iv=1;iv<=14;iv++){
        if(iv<nV){
            st1=List.at(iv).toLocal8Bit().constData();
            distCoeffs_2.at<double>(iv-1)=st1.toDouble();
        }
        else
            distCoeffs_2.at<double>(iv-1)=0.;
    }
    ui->lineEdit_fx_2->setText(QString::number(cameraMatrix_2.at<double>(0,0)));
    ui->lineEdit_fy_2->setText(QString::number(cameraMatrix_2.at<double>(1,1)));
    ui->lineEdit_cx_2->setText(QString::number(cameraMatrix_2.at<double>(0,2)));
    ui->lineEdit_cy_2->setText(QString::number(cameraMatrix_2.at<double>(1,2)));
    ui->lineEdit_k1_2->setText(QString::number(distCoeffs_2.at<double>(0)));
    ui->lineEdit_k2_2->setText(QString::number(distCoeffs_2.at<double>(1)));
    ui->lineEdit_p1_2->setText(QString::number(distCoeffs_2.at<double>(2)));
    ui->lineEdit_p2_2->setText(QString::number(distCoeffs_2.at<double>(3)));
    ui->lineEdit_k3_2->setText(QString::number(distCoeffs_2.at<double>(4)));
    ui->lineEdit_k4_2->setText(QString::number(distCoeffs_2.at<double>(5)));
    ui->lineEdit_k5_2->setText(QString::number(distCoeffs_2.at<double>(6)));
    ui->lineEdit_k6_2->setText(QString::number(distCoeffs_2.at<double>(7)));
    ui->lineEdit_s1_2->setText(QString::number(distCoeffs_2.at<double>(8)));
    ui->lineEdit_s2_2->setText(QString::number(distCoeffs_2.at<double>(9)));
    ui->lineEdit_s3_2->setText(QString::number(distCoeffs_2.at<double>(10)));
    ui->lineEdit_s4_2->setText(QString::number(distCoeffs_2.at<double>(11)));
    ui->lineEdit_tx_2->setText(QString::number(distCoeffs_2.at<double>(12)));
    ui->lineEdit_ty_2->setText(QString::number(distCoeffs_2.at<double>(13)));
    fx_2=cameraMatrix_2.at<double>(0,0);
    fy_2=cameraMatrix_2.at<double>(1,1);
    cx_2=cameraMatrix_2.at<double>(0,2);
    cy_2=cameraMatrix_2.at<double>(1,2);
    int thres;
    streamr>>comment>>thres;
    ui->sB_gL -> setValue(thres);
    streamr>>comment>>expTime;
    ui->dSBexpTime -> setValue(expTime);
    streamr>>comment>>Ps[0]>>Ps[1]>>Ps[2]>>Es[0]>>Es[1]>>Es[2];
    ui-> dSB_Xs    -> setValue(Ps[0]);
    ui-> dSB_Ys    -> setValue(Ps[1]);
    ui-> dSB_Zs    -> setValue(Ps[2]);
    ui-> dSB_Yaw_s   -> setValue(Es[0]); //yaw
    ui-> dSB_Pitch_s -> setValue(Es[1]); //pitch
    ui-> dSB_Roll_s  -> setValue(Es[2]); //roll
    for(int i=0;i<3;i++)
        Es[i]= Es[i]/180.*pig;
    streamr>>comment>>focal;
    ui-> dSB_focal -> setValue(focal);
    int iPanel;
    streamr>>comment>>iPanel;
    ui -> comboBox_panel ->setCurrentIndex(iPanel);
    double Zvalue;
    streamr>>comment>>Zvalue;
    if(iPanel==4){
        zMirror=Zvalue;
        ui->dSB_zMirror -> setValue(zMirror);
        ui->doubleSpinBox_dz0hold-> setValue(0.);
    }
    else{
        zMirror=0.;
        ui->dSB_zMirror -> setValue(0.);
        ui->doubleSpinBox_dz0hold-> setValue(Zvalue);
    }
    streamr>>comment>>Xmin>>Xmax>>Xstep;
    ui->dSB_Xmin->setValue(Xmin);
    ui->dSB_Xmax->setValue(Xmax);
    ui->dSB_Xstep->setValue(Xstep);
    streamr>>comment>>x4set;
    ui->dSB_x4set->setValue(x4set);
    XlabCam1=double(iXcur)+Pc_1[0]-x4set;
    ui->lineEdit_XlabCam1->setText(QString::number(XlabCam1));
    XlabCam2=double(iXcur)+Pc_2[0]-x4set;
    ui->lineEdit_XlabCam2->setText(QString::number(XlabCam2));

    compSou();//calcolo XYZ dei punti sorgente

    file.close();
    iReading=0;
}


void VISproPT::setWin(const std::string& _winname)
{
    cv::namedWindow(_winname);
    this->winname = _winname;
    cv::setMouseCallback(winname, on_mouse, this);
}



void VISproPT::on_mouse_internal(int event, int x, int y){
    Point point;
    //uint16_t value=0;
    //uint8_t value;
    if(x>=0 && x<Width && y>=0 && y<Height){
        ui->lineEdit_IJ -> setText(QString::number(x)+" , "+QString::number(y));
        //value=img.at<uint16_t>(y,x);
        //value=img1.at<uint8_t>(y,x);
        //ui->lineEdit_GL -> setText(QString::number(value));
    }
    if (event == EVENT_LBUTTONDOWN){
        point = Point(x, y);
        if(x<0) x=0;
        if(x>Width) x=Width;
        if(y<0) y=0;
        if(y>Height) y=Height;
        mpxX1=x;//j
        mpxY1=y;//i
        drag=1;
        printf("left button DW j=%d i=%d\n",mpxX1,mpxY1);
    }
    if (event == EVENT_LBUTTONUP){
        point = Point(x, y);
        if(x<0) x=0;
        if(x>Width) x=Width;
        if(y<0) y=0;
        if(y>Height) y=Height;
        mpxX2=x;
        mpxY2=y;
        if(mpxX2 > mpxX1 && mpxY2 > mpxY1)
            drag=2;
        else
            drag=1;
        //printf("left button UP j=%d i=%d\n",mpxX2,mpxY2);
    }
    if (event == EVENT_RBUTTONDOWN) {
        drag=-1;
        point = Point(x, y);
        if(x<0) x=0;
        if(x>Width) x=Width;
        if(y<0) y=0;
        if(y>Height) y=Height;
        mpxX1=x;
        mpxY1=y;
        //value=img.at<uint16_t>(mpxY1,mpxX1);
        //value=img.at<uint8_t>(mpxY1,mpxX1);
        printf("right button DW j=%d i=%d \n",mpxX1,mpxY1);
    }
    //    if (event == EVENT_MBUTTONDOWN) {
    //        drag=10;
    //        point = Point(x, y);
    //        if(x<0) x=0;
    //        if(x>jMAX-jMIN) x=jMAX-jMIN;
    //        if(y<0) y=0;
    //        if(y>iMAX-iMIN) y=iMAX-iMIN;
    //        mpxX1=x;
    //        mpxY1=y;
    //        value=img.at<uint16_t>(iMIN+mpxY1,jMIN+mpxX1);
    //        //value=img.at<uint8_t>(iMIN+mpxY1,jMIN+mpxX1);
    //        ui->lineEdit_IJ -> setText(QString::number(jMIN+mpxX1)+" , "+QString::number(iMIN+mpxY1));
    //        ui->lineEdit_GL -> setText(QString::number(value));
    //        printf("medium button DW j=%d i=%d grayLevel=%d\n",mpxX1,mpxY1,value);
    //    }
}


void VISproPT::selectImg(){
    filename1 = QFileDialog::getOpenFileName(this,tr("Open Image"),
                                                        dir,
                                                        tr("Image Files (* *.png *.jpg *.bmp)"));
    //cout<<"setion: "<<filename1.section('/',0,-2).toStdString()<<"\n";
    if(filename1.isEmpty())
        return;
    filename2 = QFileDialog::getOpenFileName(this,tr("Open Image"),
                                                        filename1.section('/',0,-2),
                                                        tr("Image Files (* *.png *.jpg *.bmp)"));
    if(filename2.isEmpty())
        return;
    ui->lineEdit_img -> setText(filename1+" & "+filename2);
    viewImg(filename1,filename2);
}

void VISproPT::setIntensify(){
    if(iReading==1)
        return;
    viewImg(filename1,filename2);
}


void VISproPT::setFrmCam1( ){
    setFrameCam[0][0][0]=ui->sB_jMin_1->value();
    setFrameCam[0][1][0]=ui->sB_iMin_1->value();
    setFrameCam[1][0][0]=ui->sB_wROI_1->value();
    setFrameCam[1][1][0]=ui->sB_hROI_1->value();
    QString command=part1cmd+"set_frame_cam1.cmd\t"
            +QString::number(setFrameCam[0][0][0])+"\t"
            +QString::number(setFrameCam[0][1][0])+"\t"
            +QString::number(setFrameCam[1][0][0])+"\t"
            +QString::number(setFrameCam[1][1][0]);
    system(command.toStdString().c_str());
    printf("setFrame_cam1 to %d %d %d %d\n",setFrameCam[0][0][0],setFrameCam[0][1][0],
                                            setFrameCam[1][0][0],setFrameCam[1][1][0]);
}


void VISproPT::setFrmCam2( ){
    setFrameCam[0][0][1]=ui->sB_jMin_2->value();
    setFrameCam[0][1][1]=ui->sB_iMin_2->value();
    setFrameCam[1][0][1]=ui->sB_wROI_2->value();
    setFrameCam[1][1][1]=ui->sB_hROI_2->value();
    QString command=part1cmd+"set_frame_cam2.cmd\t"
            +QString::number(setFrameCam[0][0][1])+"\t"
            +QString::number(setFrameCam[0][1][1])+"\t"
            +QString::number(setFrameCam[1][0][1])+"\t"
            +QString::number(setFrameCam[1][1][1]);
    system(command.toStdString().c_str());
    printf("setFrame_cam2 to %d %d %d %d\n",setFrameCam[0][0][1],setFrameCam[0][1][1],
                                            setFrameCam[1][0][1],setFrameCam[1][1][1]);
}


void VISproPT::grabFrame(){
    ui->lineEdit_status ->setText("grabFrame...");
    int iGrabOpt=ui->comboBox_grabOpt->currentIndex();
    QString command;
    imgN=ui->spinBox_imgN ->value();
    double newET=ui->dSBexpTime -> value();
    if(NINT(newET)!=NINT(expTime)){
        command=part1cmd+"set_expo_cam1.cmd "+QString::number(int(newET));
        system(command.toStdString().c_str());
        printf("set expTime_cam1 to %d\n",int(newET));
        command=part1cmd+"set_expo_cam2.cmd "+QString::number(int(newET));
        system(command.toStdString().c_str());
        printf("set expTime_cam2 to =%d\n",int(newET));
        expTime=newET;
    }
    if(iGrabOpt<2){
        command=part1cmd+"get_and_save_image_cam1.cmd";
        printf("command= %s\n",command.toStdString().c_str());
        system(command.toStdString().c_str());
        command="cp /mnt/data/imgcam1.bmp "+dir+"/imgCam1.bmp";
        printf("command= %s\n",command.toStdString().c_str());
        system(command.toStdString().c_str());
    }
    if(iGrabOpt!=1){
        command=part1cmd+"get_and_save_image_cam2.cmd";
        printf("command= %s\n",command.toStdString().c_str());
        system(command.toStdString().c_str());
        command="cp /mnt/data/imgcam2.bmp "+dir+"/imgCam2.bmp";
        printf("command= %s\n",command.toStdString().c_str());
        system(command.toStdString().c_str());
    }
    printf("grabbing completed\n");
    waitKey(1000);
    filename1=dir+"/imgCam1.bmp";
    filename2=dir+"/imgCam2.bmp";
    viewImg(filename1,filename2);
    Qt::CheckState state=ui-> checkBox_saveGrbImg ->checkState();
    if(state==Qt::Checked){
        QMessageBox msgBox1;
        msgBox1.setText("ATTENTION!");
        msgBox1.setInformativeText("Do you want to save the grabbed image?\n");
        msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox1.setDefaultButton(QMessageBox::Yes);
        int ret = msgBox1.exec();
        if(ret==QMessageBox::Yes){
            imgN++;
            ui->spinBox_imgN->setValue(imgN);
            if(iGrabOpt<2){
                command="cp /mnt/data/imgcam1.bmp "+dir+"/imgCam1_"+QString::number(imgN)+".bmp";
                system(command.toStdString().c_str());
            }
            if(iGrabOpt!=1){
                command="cp /mnt/data/imgcam2.bmp "+dir+"/imgCam2_"+QString::number(imgN)+".bmp";
                system(command.toStdString().c_str());
            }
        }
    }
    ui->lineEdit_status ->setText("ready");
}



void VISproPT::viewImg(QString fileImgCam1,QString fileImgCam2){
    readParam();
    printf("->viewImg(%s, %s)\n",fileImgCam1.toStdString().c_str(),fileImgCam2.toStdString().c_str());
    fflush(stdout);
    Qt::CheckState state;
    Mat img1,img2;
    Mat imgCam1=imread(fileImgCam1.toStdString(),IMREAD_ANYDEPTH);//IMREAD_UNCHANGED);
    if(! imgCam1.data){
        cout << fileImgCam1.toStdString() <<" for viewImg is not accessible!!!\n"<<"\n";
        return;
    }
    Mat imgCam2=imread(fileImgCam2.toStdString(),IMREAD_ANYDEPTH);//IMREAD_UNCHANGED);
    if(! imgCam2.data){
        cout << fileImgCam2.toStdString() <<" for viewImg is not accessible!!!\n"<<"\n";
        return;
    }
    int widthNow=imgCam1.cols;
    int heightnow=imgCam1.rows;
    if(widthNow!=Width || heightnow!=Height){
        printf("widthNow=%d <> Width=%d heightNow=%d <> Height=%d\ncall setFrmCam...\n",widthNow,Width,heightnow,Height);
        setFrmCam1();
        setFrmCam2();
        return;
    }
    state=ui-> checkBox_undist ->checkState();
    if(state==Qt::Checked)
        undistort(imgCam1,img1,cameraMatrix_1,distCoeffs_1);
    else
        img1=imgCam1;
    if(state==Qt::Checked)
        undistort(imgCam2,img2,cameraMatrix_2,distCoeffs_2);
    else
        img2=imgCam2;
    int height= img1.rows;
    int width = img1.cols;
    int thres= ui->sB_gL -> value();   //threshold
    Mat matDisplay8bCam1(height,width,CV_8UC3);
    Mat matDisplay8bCam2(height,width,CV_8UC3);
    uint8_t value8cam1,value8cam2;
    for(int i=0; i<height; i++){
        for(int j=0; j< width; j++){
            value8cam1=img1.at<uint8_t>(i,j);//8 bit
            value8cam2=img2.at<uint8_t>(i,j);//8 bit
            for(int k=0; k<3; k++){
                if(value8cam1<thres)
                    matDisplay8bCam1.at<Vec3b>(i,j)[k]=value8cam1;
                else{
                    if(k<2)
                       matDisplay8bCam1.at<Vec3b>(i,j)[k]=0;
                    else
                       matDisplay8bCam1.at<Vec3b>(i,j)[k]=255;
                }
                if(value8cam2<thres)
                    matDisplay8bCam2.at<Vec3b>(i,j)[k]=value8cam2;
                else{
                    if(k<2)
                       matDisplay8bCam2.at<Vec3b>(i,j)[k]=0;
                    else
                       matDisplay8bCam2.at<Vec3b>(i,j)[k]=255;
                }
            }
        }
    }
    imshow("ROIcam1",img1);
    //mpxX1=x;//j
    //mpxY1=y;//i
    //rectangle(matDisplay8bCam1, Pt1,Pt2,Scalar(255),1, 8, 0 );
    Pt1.x=mpxX1;
    Pt1.y=mpxY1;
    if(mpxX1>=0 && mpxX1<Width && mpxY1>=0 && mpxY1<Height){
        circle(matDisplay8bCam1, Pt1,2, Vec3b(0,255,0), -1);
        imshow("Camera1",matDisplay8bCam1);
        unitV12(0.0,0.0,0.0,(mpxX1-cx_1)/fx_1,(mpxY1-cy_1)/fy_1,1.);
        Plane2Earth(EA_1[0],EA_1[1],EA_1[2],Vservice[0],Vservice[1],Vservice[2]);
        double dZcam=0.,dPitch=0.,dRoll=0.;
        state=ui->checkBox_artCor->checkState();
        if(state==Qt::Checked){
            int iLastArtCor=215;
            while(corArte[iLastArtCor][0]>iXcur && iLastArtCor>0){
                iLastArtCor--;
            }
            dZcam=corArte[iLastArtCor][1];
            dPitch=corArte[iLastArtCor][3];
            dRoll=corArte[iLastArtCor][4];
            printf("xNear=%f dZcam=%f dPitch=%f dRoll=%f\n",
                   corArte[iLastArtCor][0],dZcam,dPitch,dRoll);

        }
        Plane2Earth(0.,-dPitch,-dRoll,Vservice[0],Vservice[1],Vservice[2]);
        double disto1=abs((Pc_1[2]+dZcam)/Vservice[2]);
        printf("Cam1 (x,y,z)_dot= %f %f %f\n",
               double(iXcur)+Pc_1[0]-x4set+disto1*Vservice[0],
                Pc_1[1]+disto1*Vservice[1],
                Pc_1[2]+dZcam+disto1*Vservice[2]);
    }
    imshow("ROIcam2",img2);
    //rectangle(matDisplay8bCam2, Pt1,Pt2,Scalar(255),1, 8, 0 );
    Pt1.x=cx_2;
    Pt1.y=cy_2;
    circle(matDisplay8bCam2, Pt1, 10, Vec3b(0,255,0), -1);
    imshow("Camera2",matDisplay8bCam2);
    Plane2Earth(EA_2[0],EA_2[1],EA_2[2],0.,0.,1.);
    double disto2=abs(Pc_2[2]/Vservice[2]);
    printf("Cam2 (x,y,z)_dot= %f %f %f\n",
           Pc_2[0]+disto2*Vservice[0],
           Pc_2[1]+disto2*Vservice[1],
           Pc_2[2]+disto2*Vservice[2]);


    state=ui-> checkBox_make8bRGBcopy ->checkState();
    if(state==Qt::Checked){
        imwrite(dir.toStdString()+"/matDisplay8bCam1.jpg",matDisplay8bCam1);
        imwrite(dir.toStdString()+"/matDisplay8bCam2.jpg",matDisplay8bCam2);
    }

    waitKey(10);
    imgCam1.release();
    imgCam2.release();
    img1.release();
    img2.release();
    matDisplay8bCam1.release();
    matDisplay8bCam2.release();
}


void VISproPT::setCB(){
    int iChessType=ui->comboBox_chessType->currentIndex();
    if(iChessType==0){
        ui-> sB_nCornHor  -> setValue(37);
        ui-> sB_nCornVer  -> setValue(17);
        ui->dSB_cornerStep-> setValue(24.9);
    }
    else if(iChessType==1){
        ui-> sB_nCornHor  -> setValue(15);
        ui-> sB_nCornVer  -> setValue(9);
        ui->dSB_cornerStep-> setValue(69.83);
    }
    else if(iChessType==2){
        ui-> sB_nCornHor  -> setValue(26);
        ui-> sB_nCornVer  -> setValue(21);
        ui->dSB_cornerStep-> setValue(19.94);
    }
    else if(iChessType==3){
        ui-> sB_nCornHor  -> setValue(53);
        ui-> sB_nCornVer  -> setValue(45);
        ui->dSB_cornerStep-> setValue(10.01);
    }
    else{
        ui-> sB_nCornHor  -> setValue(45);
        ui-> sB_nCornVer  -> setValue(29);
        ui->dSB_cornerStep-> setValue(24.98);
    }
}



void VISproPT::CALLcalib1(){
    calibrate(1);
}

void VISproPT::CALLcalib2(){
    calibrate(2);
}


void VISproPT::calibrate(int iCam){
    readParam();
    ui->lineEdit_status ->setText("calibrate Cam"+QString::number(iCam));
    Qt::CheckState state;
    state=ui-> checkBox_make8bRGBcopy ->checkState();

    int numCornersHor = ui-> sB_nCornHor -> value();
    int numCornersVer = ui-> sB_nCornVer -> value();
    int numSquares = numCornersHor * numCornersVer;
    double cornerStep = ui->dSB_cornerStep -> value();
    Size board_sz = Size(numCornersHor, numCornersVer);
    vector<vector<Point3f>> object_points;
    vector<vector<Point2f>> image_points;
    vector<Point2f> corners;
    Mat imageCal8b;
    Mat imageCal(Height,Width,CV_8UC1);
    Mat imageRGB(Height,Width,CV_8UC3);
    uint8_t val8;
    vector<Point3f> obj;
    for(int j=0;j<numSquares;j++){
        //obj.push_back(Point3d((j/numCornersHor)*cornerStep, (j%numCornersHor)*cornerStep, 0.0));
        obj.push_back(Point3d((j%numCornersHor)*cornerStep,(j/numCornersHor)*cornerStep,0.0));
        //printf("obj[%d]=(%f,%f,%f)\n",j,obj[j].x,obj[j].y,obj[j].z);
    }
//    for(int i=0;i<numCornersVer;i++){
//        for(int j=0;j<numCornersHor;j++)
//          printf("Point[%d][%d]=(%f,%f,%f)\n",j,i,obj[i*numCornersVer+j].x,
//                  obj[i*numCornersVer+j].y,obj[i*numCornersVer+j].z);
//    }
    QString fileCalImg,file2save;
    printf("Camera Calibration is started...\n");
    PathPro="";
    PathPro=QFileDialog::getExistingDirectory(this, tr("Open Directory for calibrate"),
                                              dir+"/camCalib",
                                              QFileDialog::ShowDirsOnly);
    if(PathPro.isEmpty()){
        printf("abort!\n");
        return;
    }
    PathPro=PathPro+"/";
    QMessageBox msgBox1;
    int iOK=0;
    int iContinue=1;
    int suc=1;
    int sucLim=ui->spinBox_imgN->value();
    while(iContinue==1){
        fileCalImg=PathPro+"imgCam"+QString::number(iCam)+"_"+QString::number(suc)+".bmp";
        cout << "processing img "<<fileCalImg.toStdString() <<"\n";
        QFileInfo check_file(fileCalImg);
        // check if file exists and if yes: Is it really a file and no directory?
        if (check_file.exists() && check_file.isFile()) {
            imageCal8b=imread(fileCalImg.toStdString(),IMREAD_ANYDEPTH);
        } else {
            iContinue=0;
        }
        if(iContinue==1){
            for(int i=0; i<Height; i++){
                for(int j=0; j< Width; j++){
                    val8=imageCal8b.at<uint8_t>(i,j);
                    imageCal.at<uint8_t>(i,j)=val8;
                    for(int k=0; k<3; k++)
                        imageRGB.at<Vec3b>(i,j)[k]=val8;
                }
            }
            imshow("Board",imageRGB);
            if(state==Qt::Checked){
                file2save=PathPro+"img_current"+QString::number(suc)+".jpg";
                imwrite(file2save.toStdString(),imageCal);
            }
            waitKey(10);
            bool found = findChessboardCorners(imageCal, board_sz, corners, CALIB_CB_ADAPTIVE_THRESH | CALIB_CB_FILTER_QUADS);
            if(found){
                cornerSubPix(imageCal, corners, Size(11, 11), Size(-1, -1), TermCriteria( TermCriteria::EPS+TermCriteria::COUNT, 30, 0.0001 ));
                drawChessboardCorners(imageRGB, board_sz, corners, found);
            }
            imshow("Board",imageRGB);
            msgBox1.setText("ATTENTION!");
            msgBox1.setInformativeText("Were corners rightly located?\n");
            msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox::Abort);
            msgBox1.setDefaultButton(QMessageBox::Yes);
            int ret = msgBox1.exec();
            if(ret==QMessageBox::Yes){
                image_points.push_back(corners);
                object_points.push_back(obj);
                cout <<"positive!\n";
//                for(int iCor=0;iCor<numSquares;iCor++)
//                    printf("%d %f %f\n",iCor,corners[iCor].x,corners[iCor].y);
                iOK++;
            }
            else if(ret==QMessageBox::No)
                cout <<"negative!\n";
            else if(ret==QMessageBox::Abort)
                iContinue=0;
            suc++;
            printf("iOK= %d\n",iOK);
        }

        if(suc==sucLim)iContinue=0;
    }
    if(iOK>0){
        cout<<"calibreCamera...."<<"\n";
        vector<Mat> rvecs;
        vector<Mat> tvecs;
        Mat corrImg;
        int iMth=ui->comboBox_calibMethod->currentIndex();
        int method=0;
        int Ndim=5;
        if(iMth==0){
            method=CALIB_FIX_ASPECT_RATIO;
            Ndim=5;
        }
        else if(iMth==1){
            method=CALIB_FIX_ASPECT_RATIO+CALIB_ZERO_TANGENT_DIST+
                    CALIB_FIX_K1+CALIB_FIX_K2+CALIB_FIX_K3;
                    //CALIB_FIX_K4+CALIB_FIX_K5+CALIB_FIX_K6+
                    //CALIB_FIX_S1_S2_S3_S4+
                    //CALIB_FIX_TAUX_TAUY;
            Ndim=5;
        }
        else if(iMth==2){
            method=CALIB_FIX_ASPECT_RATIO+CALIB_RATIONAL_MODEL;
            Ndim=8;
        }
        else if(iMth==3){
            method=CALIB_FIX_ASPECT_RATIO+CALIB_THIN_PRISM_MODEL;
            Ndim=12;
        }
        else if(iMth==4){
            method=CALIB_FIX_ASPECT_RATIO+CALIB_TILTED_MODEL;
            Ndim=14;
        }
        else if(iMth==5){
            method=CALIB_FIX_ASPECT_RATIO+CALIB_RATIONAL_MODEL+
                    CALIB_THIN_PRISM_MODEL+CALIB_TILTED_MODEL;
            Ndim=14;
        }
        else if(iMth==6){
            method=CALIB_USE_INTRINSIC_GUESS+CALIB_FIX_FOCAL_LENGTH+
                    CALIB_FIX_PRINCIPAL_POINT;
            Ndim=5;
        }
        else{
            method=CALIB_USE_INTRINSIC_GUESS+CALIB_FIX_FOCAL_LENGTH+
                    CALIB_FIX_PRINCIPAL_POINT+CALIB_RATIONAL_MODEL+
                    CALIB_THIN_PRISM_MODEL+CALIB_TILTED_MODEL;
            Ndim=14;
        }
        printf("method=%d\n",method);

        if(iCam==1){
            Mat discoe;
            //if(iMth<=5)
                calibrateCamera(object_points, image_points, imageCal.size(), cameraMatrix_1, discoe, rvecs, tvecs,method);
            //else
            //    solvePnP(object_points, image_points, cameraMatrix_1, discoe, rvecs, tvecs,false,SOLVEPNP_ITERATIVE);
            for(int j=0;j<Ndim;j++){
                distCoeffs_1.at<double>(j)=discoe.at<double>(j);
            }
            for(int j=Ndim;j<14;j++){
                distCoeffs_1.at<double>(j)=0.;
            }
            for(int j=0;j<14;j++){
                if(abs(distCoeffs_1.at<double>(j))<1.e-100)
                    distCoeffs_1.at<double>(j)=0.;
            }
            ui->lineEdit_fx_1->setText(QString::number(cameraMatrix_1.at<double>(0,0)));
            ui->lineEdit_fy_1->setText(QString::number(cameraMatrix_1.at<double>(1,1)));
            ui->lineEdit_cx_1->setText(QString::number(cameraMatrix_1.at<double>(0,2)));
            ui->lineEdit_cy_1->setText(QString::number(cameraMatrix_1.at<double>(1,2)));
            ui->lineEdit_k1_1->setText(QString::number(distCoeffs_1.at<double>(0)));
            ui->lineEdit_k2_1->setText(QString::number(distCoeffs_1.at<double>(1)));
            ui->lineEdit_p1_1->setText(QString::number(distCoeffs_1.at<double>(2)));
            ui->lineEdit_p2_1->setText(QString::number(distCoeffs_1.at<double>(3)));
            ui->lineEdit_k3_1->setText(QString::number(distCoeffs_1.at<double>(4)));
            ui->lineEdit_k4_1->setText(QString::number(distCoeffs_1.at<double>(5)));
            ui->lineEdit_k5_1->setText(QString::number(distCoeffs_1.at<double>(6)));
            ui->lineEdit_k6_1->setText(QString::number(distCoeffs_1.at<double>(7)));
            ui->lineEdit_s1_1->setText(QString::number(distCoeffs_1.at<double>(8)));
            ui->lineEdit_s2_1->setText(QString::number(distCoeffs_1.at<double>(9)));
            ui->lineEdit_s3_1->setText(QString::number(distCoeffs_1.at<double>(10)));
            ui->lineEdit_s4_1->setText(QString::number(distCoeffs_1.at<double>(11)));
            ui->lineEdit_tx_1->setText(QString::number(distCoeffs_1.at<double>(12)));
            ui->lineEdit_ty_1->setText(QString::number(distCoeffs_1.at<double>(13)));
            fx_1=cameraMatrix_1.at<double>(0,0);
            fy_1=cameraMatrix_1.at<double>(1,1);
            cx_1=cameraMatrix_1.at<double>(0,2);
            cy_1=cameraMatrix_1.at<double>(1,2);
            fcam_1=(cameraMatrix_1.at<double>(0,0)*pxdimX+cameraMatrix_1.at<double>(1,1)*pxdimY)/2.;
            ui-> dSB_fcam_1  -> setValue(fcam_1);
            undistort(imageRGB,corrImg,cameraMatrix_1,distCoeffs_1);
            QFile fLog(dir+"/fLog.txt");
            if (!fLog.open(QIODevice::WriteOnly | QIODevice::Append))
                return;
            QTextStream out(&fLog);
            out<<iCam<<"\t"<<suc-1<<"\t"<<cameraMatrix_1.at<double>(0,0)<<"\t"<<cameraMatrix_1.at<double>(1,1)<<"\t"
              <<cameraMatrix_1.at<double>(0,2)<<"\t"<<cameraMatrix_1.at<double>(1,2)<<"\n";
            fLog.close();
        }
        else{
            Mat discoe;
            //if(iMth<=5)
                calibrateCamera(object_points, image_points, imageCal.size(), cameraMatrix_2, discoe, rvecs, tvecs,method);
            //else
            //    solvePnP(object_points, image_points, cameraMatrix_2, discoe, rvecs, tvecs,false,SOLVEPNP_ITERATIVE);
            for(int j=0;j<Ndim;j++){
                distCoeffs_2.at<double>(j)=discoe.at<double>(j);
            }

            for(int j=Ndim;j<14;j++){
                distCoeffs_2.at<double>(j)=0.;
            }
            for(int j=0;j<14;j++){
                if(abs(distCoeffs_2.at<double>(j))<1.e-100)
                    distCoeffs_2.at<double>(j)=0.;
            }
            ui->lineEdit_fx_2->setText(QString::number(cameraMatrix_2.at<double>(0,0)));
            ui->lineEdit_fy_2->setText(QString::number(cameraMatrix_2.at<double>(1,1)));
            ui->lineEdit_cx_2->setText(QString::number(cameraMatrix_2.at<double>(0,2)));
            ui->lineEdit_cy_2->setText(QString::number(cameraMatrix_2.at<double>(1,2)));
            ui->lineEdit_k1_2->setText(QString::number(distCoeffs_2.at<double>(0)));
            ui->lineEdit_k2_2->setText(QString::number(distCoeffs_2.at<double>(1)));
            ui->lineEdit_p1_2->setText(QString::number(distCoeffs_2.at<double>(2)));
            ui->lineEdit_p2_2->setText(QString::number(distCoeffs_2.at<double>(3)));
            ui->lineEdit_k3_2->setText(QString::number(distCoeffs_2.at<double>(4)));
            ui->lineEdit_k4_2->setText(QString::number(distCoeffs_2.at<double>(5)));
            ui->lineEdit_k5_2->setText(QString::number(distCoeffs_2.at<double>(6)));
            ui->lineEdit_k6_2->setText(QString::number(distCoeffs_2.at<double>(7)));
            ui->lineEdit_s1_2->setText(QString::number(distCoeffs_2.at<double>(8)));
            ui->lineEdit_s2_2->setText(QString::number(distCoeffs_2.at<double>(9)));
            ui->lineEdit_s3_2->setText(QString::number(distCoeffs_2.at<double>(10)));
            ui->lineEdit_s4_2->setText(QString::number(distCoeffs_2.at<double>(11)));
            ui->lineEdit_tx_2->setText(QString::number(distCoeffs_2.at<double>(12)));
            ui->lineEdit_ty_2->setText(QString::number(distCoeffs_2.at<double>(13)));
            fx_2=cameraMatrix_2.at<double>(0,0);
            fy_2=cameraMatrix_2.at<double>(1,1);
            cx_2=cameraMatrix_2.at<double>(0,2);
            cy_2=cameraMatrix_2.at<double>(1,2);
            fcam_2=(cameraMatrix_2.at<double>(0,0)*pxdimX+cameraMatrix_2.at<double>(1,1)*pxdimY)/2.;
            ui-> dSB_fcam_2  -> setValue(fcam_2);
            undistort(imageRGB,corrImg,cameraMatrix_2,distCoeffs_2);
            QFile fLog(dir+"/fLog.txt");
            if (!fLog.open(QIODevice::WriteOnly | QIODevice::Append))
                return;
            QTextStream out(&fLog);
            out<<iCam<<"\t"<<suc-1<<"\t"<<cameraMatrix_2.at<double>(0,0)<<"\t"<<cameraMatrix_2.at<double>(1,1)<<"\t"
              <<cameraMatrix_2.at<double>(0,2)<<"\t"<<cameraMatrix_2.at<double>(1,2)<<"\n";
            fLog.close();
        }
        msgBox1.setText("ATTENTION!");
        msgBox1.setInformativeText("Do you want to save the new values?\n");
        msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox1.setDefaultButton(QMessageBox::Yes);
        int ret = msgBox1.exec();
        if(ret==QMessageBox::Yes)
            callSaveAcam();
        else
            readAcam(fileAcam);
        cout<<"done!\n"<<"\n";
    }
    else{
        msgBox1.setText("ATTENTION!");
        msgBox1.setInformativeText("Number of succesful board is 0\n");
        msgBox1.setStandardButtons(QMessageBox::Yes);
        msgBox1.setDefaultButton(QMessageBox::Yes);
        msgBox1.exec();
    }
    imageCal.release();
    imageCal8b.release();
    imageRGB.release();
    ui->lineEdit_status ->setText("ready");
    fflush(stdout);
}


void VISproPT::refBFM(){
    Qt::CheckState state;
    state=ui->checkBox_artCor->checkState();
    iArtCor=0;
    if(state==Qt::Checked)
        iArtCor=1;
    int nFrame;
    double value;
    if(iTot<1)
    dirMeas = QFileDialog::getExistingDirectory(this, tr("Open Directory for process"),
                                                        dir,
                                                        QFileDialog::ShowDirsOnly);
    printf("dirMeas= %s\n",dirMeas.toStdString().c_str());
    Xmax=ui->dSB_Xmax->value();
    Xmin=ui->dSB_Xmin->value();
    Xstep=ui->dSB_Xstep->value();
    Xstart=Xmax;
    nFrame=NINT((Xmax-Xmin)/Xstep)+1;
    Xstep=-Xstep;
    Dx21=ui->dSB_Dx21->value();
    Dz21=ui->dSB_Dz21->value();
    printf("nFrame_expected=%d xStep=%f\n",nFrame,Xstep);
    int iRefOpt=ui->comboBox_refOpt->currentIndex();

    // set up the parameters (check the defaults in opencv's code in blobdetector.cpp)
    //cv::SimpleBlobDetector::Params params;

    // Filter by treshold
    params.thresholdStep=5;
    params.minThreshold= ui -> sB_gL -> value();
    params.maxThreshold=255;

    //Filter by distance
    minDis=ui->doubleSpinBox_minDist->value();
    params.minDistBetweenBlobs = float(minDis);

    // Filter by Area
    state=ui->checkBox_Area->checkState();
    if(state==Qt::Checked){
        params.filterByArea = true;
        value=ui->doubleSpinBox_minArea->value();
        params.minArea = float(value);
        params.maxArea = 500.0f;
    }
    else
        params.filterByArea = false;

    // Filter by Circularity
    state=ui->checkBox_circularity->checkState();
    if(state==Qt::Checked){
        params.filterByCircularity = true;
        value=ui->doubleSpinBox_minCircularity->value();
        params.minCircularity =float(value);
    }
    else
        params.filterByCircularity = false;

    //filter by Convexity
    state=ui->checkBox_convexity->checkState();
    if(state==Qt::Checked){
        params.filterByConvexity = true;
        value=ui->doubleSpinBox_minConvexity->value();
        params.minConvexity =float(value);
    }
    else
        params.filterByConvexity = false;

    // Filter by Inertia
    state=ui->checkBox_inertia->checkState();
    if(state==Qt::Checked){
        params.filterByInertia = true;
        value=ui->doubleSpinBox_minInertia->value();
        params.minInertiaRatio =float(value);
    }
    else
        params.filterByInertia = false;

    params.filterByColor = false;
    // Filter by Circularity
    params.filterByCircularity = false;
    //params.minCircularity = 0.8f;
    int info;
    int m=2*(isMAX-isMIN-9)*(jsMAX-jsMIN-9);
    int n;
    if(iRefOpt==0)
        n=6;
    else if(iRefOpt==1)
        n=10;
    else
        n=16;
    int lwa=m*n+5*n+m;
    int* iwa=nullptr;
    iwa = new int[n];
    double* xx=nullptr;
    xx = new double[n];
    double* fvec=nullptr;
    fvec = new double[m];
    double* wa=nullptr;
    wa = new double[lwa];
    double tol=sqrt(dpmpar(1));
    tol=10.*tol;
    pTF2[0].Nt=iRefOpt;

    if(iRefOpt==0){
        xx[0]=fx_1;
        xx[1]=cx_1;
        xx[2]=cy_1;
        xx[3]=fx_2;
        xx[4]=cx_2;
        xx[5]=cy_2;
    }
    else if(iRefOpt==1){
        xx[0]=Pc_1[0];
        xx[1]=Pc_1[1];
        xx[2]=Pc_1[2];
        xx[3]=EA_1[0];
        xx[4]=EA_1[1];
        xx[5]=EA_1[2];
        //xx[6]=Pc_2[0];
        Pc_2[0]=Pc_1[0]+Dx21;
        xx[6]=Pc_2[1];
        //xx[8]=Pc_2[2];
        Pc_2[2]=Pc_1[2]+Dz21;
        xx[7]=EA_2[0];
        xx[8]=EA_2[1];
        xx[9]=EA_2[2];
    }
    else{
        xx[0]=fx_1;
        xx[1]=cx_1;
        xx[2]=cy_1;
        xx[3]=fx_2;
        xx[4]=cx_2;
        xx[5]=cy_2;
        xx[6]=Pc_1[0];
        xx[7]=Pc_1[1];
        xx[8]=Pc_1[2];
        xx[9]=EA_1[0];
        xx[10]=EA_1[1];
        xx[11]=EA_1[2];
        //xx[12]=Pc_2[0];
        Pc_2[0]=Pc_1[0]+Dx21;
        xx[12]=Pc_2[1];
        //xx[14]=Pc_2[2];
        Pc_2[2]=Pc_1[2]+Dz21;
        xx[13]=EA_2[0];
        xx[14]=EA_2[1];
        xx[15]=EA_2[2];
    }

    printf("refCam1&2: n=%d parameters and m=%d data\n",n,m);
    if(iRefOpt==0 || iRefOpt==2)
        printf("refCam initial values f_1=%f cx_1=%f cy_1=%f f_2=%f cx_2=%f cy_2=%f\n",
           xx[0],xx[1],xx[2],xx[3],xx[4],xx[5]);
    if(iRefOpt>0){
        printf("\tPcam1= %f\t%f\t%f\tEcam1= %f\t%f\t%f\n",
               Pc_1[0],Pc_1[1],Pc_1[2],EA_1[0]*180./pig,EA_1[1]*180./pig,EA_1[2]*180./pig);
        printf("\tPcam2= %f\t%f\t%f\tEcam1= %f\t%f\t%f\n",
               Pc_2[0],Pc_2[1],Pc_2[2],EA_2[0]*180./pig,EA_2[1]*180./pig,EA_2[2]*180./pig);
    }

    if(iTot<1){
        iTot=-1;      //results[] initialization
        kernelpro(1); //new image processing
        kernelpro(2);
    }

    info=lmdif1(fcnRefCam, &pTF2, m, n, xx,fvec, tol, iwa, wa, lwa);

    if(iRefOpt==0){
        fx_1=xx[0];
        fy_1=xx[0];
        cx_1=xx[1];
        cy_1=xx[2];
        fx_2=xx[3];
        fy_2=xx[3];
        cx_2=xx[4];
        cy_2=xx[5];
    }
    else if(iRefOpt==1){
        Pc_1[0]=xx[0];
        Pc_1[1]=xx[1];
        Pc_1[2]=xx[2];
        EA_1[0]=xx[3];
        EA_1[1]=xx[4];
        EA_1[2]=xx[5];
        //Pc_2[0]=xx[6];
        Pc_2[0]=Pc_1[0]+Dx21;
        Pc_2[1]=xx[6];
        //Pc_2[2]=xx[8];
        Pc_2[2]=Pc_1[2]+Dz21;
        EA_2[0]=xx[7];
        EA_2[1]=xx[8];
        EA_2[2]=xx[9];
    }
    else{
        fx_1=xx[0];
        fy_1=xx[0];
        cx_1=xx[1];
        cy_1=xx[2];
        fx_2=xx[3];
        fy_2=xx[3];
        cx_2=xx[4];
        cy_2=xx[5];
        Pc_1[0]=xx[6];
        Pc_1[1]=xx[7];
        Pc_1[2]=xx[8];
        EA_1[0]=xx[9];
        EA_1[1]=xx[10];
        EA_1[2]=xx[11];
        //Pc_2[0]=xx[12];
        Pc_2[0]=Pc_1[0]+Dx21;
        Pc_2[1]=xx[12];
        //Pc_2[2]=xx[14];
        Pc_2[2]=Pc_1[2]+Dz21;
        EA_2[0]=xx[13];
        EA_2[1]=xx[14];
        EA_2[2]=xx[15];
    }
    cameraMatrix_1.at<double>(0,0)=fx_1;
    cameraMatrix_1.at<double>(1,1)=fy_1;
    cameraMatrix_1.at<double>(0,2)=cx_1;
    cameraMatrix_1.at<double>(1,2)=cy_1;
    ui->lineEdit_fx_1->setText(QString::number(cameraMatrix_1.at<double>(0,0)));
    ui->lineEdit_fy_1->setText(QString::number(cameraMatrix_1.at<double>(1,1)));
    ui->lineEdit_cx_1->setText(QString::number(cameraMatrix_1.at<double>(0,2)));
    ui->lineEdit_cy_1->setText(QString::number(cameraMatrix_1.at<double>(1,2)));
    ui-> dSB_Xc_1 -> setValue(Pc_1[0]);
    ui-> dSB_Yc_1 -> setValue(Pc_1[1]);
    ui-> dSB_Zc_1 -> setValue(Pc_1[2]);
    ui-> dSB_Yaw_1   -> setValue(EA_1[0]*180./pig); //yaw
    ui-> dSB_Pitch_1 -> setValue(EA_1[1]*180./pig); //pitch
    ui-> dSB_Roll_1  -> setValue(EA_1[2]*180./pig); //roll
    XlabCam1=Pc_1[0];
    ui->lineEdit_XlabCam1->setText(QString::number(XlabCam1));
    cameraMatrix_2.at<double>(0,0)=fx_2;
    cameraMatrix_2.at<double>(1,1)=fy_2;
    cameraMatrix_2.at<double>(0,2)=cx_2;
    cameraMatrix_2.at<double>(1,2)=cy_2;
    ui->lineEdit_fx_2->setText(QString::number(cameraMatrix_2.at<double>(0,0)));
    ui->lineEdit_fy_2->setText(QString::number(cameraMatrix_2.at<double>(1,1)));
    ui->lineEdit_cx_2->setText(QString::number(cameraMatrix_2.at<double>(0,2)));
    ui->lineEdit_cy_2->setText(QString::number(cameraMatrix_2.at<double>(1,2)));
    ui-> dSB_Xc_2 -> setValue(Pc_2[0]);
    ui-> dSB_Yc_2 -> setValue(Pc_2[1]);
    ui-> dSB_Zc_2 -> setValue(Pc_2[2]);
    ui-> dSB_Yaw_2   -> setValue(EA_2[0]*180./pig); //yaw
    ui-> dSB_Pitch_2 -> setValue(EA_2[1]*180./pig); //pitch
    ui-> dSB_Roll_2  -> setValue(EA_2[2]*180./pig); //roll
    XlabCam2=Pc_2[0];
    ui->lineEdit_XlabCam2->setText(QString::number(XlabCam2));
    ui->dSB_Xmax->setValue(Xstart);
    ui->dSB_Xstep->setValue(abs(Xstep));

    if(iRefOpt==0)
        printf("refCam terminated with values f_1=%f cx_1=%f cy_1=%f f_2=%f cx_2=%f cy_2=%f\n",
           xx[0],xx[1],xx[2],xx[3],xx[4],xx[5]);
    else{
        printf("refCam terminated:\n\tPcam1= %f\t%f\t%f\tEcam1= %f\t%f\t%f\n",
               Pc_1[0],Pc_1[1],Pc_1[2],EA_1[0]*180./pig,EA_1[1]*180./pig,EA_1[2]*180./pig);
        printf("\tPcam2= %f\t%f\t%f\tEcam1= %f\t%f\t%f\n",
               Pc_2[0],Pc_2[1],Pc_2[2],EA_2[0]*180./pig,EA_2[1]*180./pig,EA_2[2]*180./pig);
    }
    ui->lineEdit_chi2->setText(QString::number(chi2fin));
    QMessageBox msgBox1;
    msgBox1.setText("ATTENTION!");
    msgBox1.setInformativeText("Do you want to save the new values?\n");
    msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox1.setDefaultButton(QMessageBox::Yes);
    int ret = msgBox1.exec();
    if(ret==QMessageBox::Yes)
        callSaveAcam();
    else
        readAcam(fileAcam);

    if(wa){
        delete[] wa;
        wa = nullptr;
    }
    if(fvec){
        delete [] fvec;
        fvec=nullptr;
    }
    if(xx){
        delete [] xx;
        xx=nullptr;
    }
    if(iwa){
        delete []  iwa;
        iwa=nullptr;
    }
    printf("refCam exit!\n");
    fflush(stdout);
}


void VISproPT::GoTo(){
    QString command;
    double X=ui->dSB_X ->value();
    ui->lineEdit_status ->setText("goto x= "+QString::number(X));
    command=part1cmd+"go_to_position_mm.cmd "+QString::number(X);
    system(command.toStdString().c_str());
    do{
        waitKey(1000);
        iXcur=getPosition();
        ui->lineEdit_Xcur->setText(QString::number(iXcur));
        XlabCam1=double(iXcur)+Pc_1[0]-x4set;
        ui->lineEdit_XlabCam1->setText(QString::number(XlabCam1));
        XlabCam2=double(iXcur)+Pc_2[0]-x4set;
        ui->lineEdit_XlabCam2->setText(QString::number(XlabCam2));
    }while(abs(X-iXcur)>=1.);
    ui->lineEdit_status ->setText("ready");
}


void VISproPT::setCamera(){
    QMessageBox msgBox1,msgBox;
    msgBox1.setText("ATTENTION!");
    msgBox1.setInformativeText("Before setting cameras you must have already load the right setting!\nIt is the case?");
    msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox1.setDefaultButton(QMessageBox::Yes);
    int ret = msgBox1.exec();
    if(ret==QMessageBox::No)
        return;
    //readParam();
    ui->lineEdit_status ->setText("setCamera");
    int iSetCamMeth=ui->comboBox_setCamMeth->currentIndex();
    PathPro="";
    PathPro=QFileDialog::getExistingDirectory(this, tr("Open Directory 4 setCamera"),
                                              dir+"/SetCam",
                                              QFileDialog::ShowDirsOnly);
    if(PathPro.isEmpty()){
        printf("abort!\n");
        return;
    }
    PathPro=PathPro+"/";
    Qt::CheckState state;
    Mat img1,img2,imageCal8b;
    Mat imgCam1=imread((PathPro+"imgCam1.bmp").toStdString(),IMREAD_ANYDEPTH);//IMREAD_UNCHANGED);
    if(!imgCam1.data){
        cout << (PathPro+"imgCam1.bmp").toStdString() <<" for viewImg is not accessible!!!\n"<<"\n";
        return;
    }
    Mat imgCam2=imread((PathPro+"imgCam2.bmp").toStdString(),IMREAD_ANYDEPTH);//IMREAD_UNCHANGED);
    if(!imgCam2.data){
        cout << (PathPro+"imgCam2.bmp").toStdString() <<" for viewImg is not accessible!!!\n"<<"\n";
        return;
    }
    state=ui-> checkBox_undist ->checkState();
    if(state==Qt::Checked)
        undistort(imgCam1,img1,cameraMatrix_1,distCoeffs_1);
    else
        img1=imgCam1;
    if(state==Qt::Checked)
        undistort(imgCam2,img2,cameraMatrix_2,distCoeffs_2);
    else
        img2=imgCam2;

    int numCornersHor = ui-> sB_nCornHor -> value();
    int numCornersVer = ui-> sB_nCornVer -> value();
    int numSquares = numCornersHor * numCornersVer;
    Size board_sz = Size(numCornersHor, numCornersVer);
    double cornerStep = ui-> dSB_cornerStep -> value();
    vector<Point2f> corners;
    double xOffset=0.;//-0.24;
    double yOffset=0.;//-1.24;
    int Nt=numSquares;
    Ncorn=Nt;
    int Np=0;
    for(int iCam=1;iCam<3;iCam++){
        if(iCam==1)
            Np=0;
        else
            Np=Nt;
        int iCor=0;
        for(int i=0; i<numCornersVer; i++){
            for(int j=0; j<numCornersHor; j++){
                P[Np+iCor][0]=cornerStep*((numCornersVer+1)/2-1-i)+xOffset;
                P[Np+iCor][1]=cornerStep*(((numCornersHor+1)/2-1)-j)+yOffset;
                P[Np+iCor][2]=0.;//z=0 is set on the top side of the 6mm thick reflective float glass
                //printf("P[%d]=Pcorn[%d][%d]= %f %f %f\n",Np+iCor,j,i,P[Np+iCor][0],P[Np+iCor][1],P[Np+iCor][2]);
                iCor++;
            }
        }
    }
    for(int iCam=1;iCam<3;iCam++){
        if(iCam==1){
            imageCal8b=img1;
            Np=0;
        }
        else{
            imageCal8b=img2;
            Np=Nt;
        }
        Mat imageRGB(Height,Width,CV_8UC3);
        uint8_t val8;
        for(int i=0; i<Height; i++){
            for(int j=0; j< Width; j++){
                val8=imageCal8b.at<uint8_t>(i,j);
                for(int k=0; k<3; k++)
                    imageRGB.at<Vec3b>(i,j)[k]=val8;
            }
        }
        imshow("Board",imageRGB);
        waitKey(10);
        cout <<"automatic search of corners ....\n" ;
        bool found = findChessboardCorners(imageCal8b, board_sz, corners, CALIB_CB_ADAPTIVE_THRESH | CALIB_CB_FILTER_QUADS);
        if(found){
            cornerSubPix(imageCal8b, corners, Size(11, 11), Size(-1, -1), TermCriteria( TermCriteria::EPS+TermCriteria::COUNT, 30, 0.0001 ));
            drawChessboardCorners(imageRGB, board_sz, corners, found);
        }
        waitKey(10);
        imshow("Board",imageRGB);
        img=imageCal8b;
        msgBox1.setText("ATTENTION!");
        msgBox1.setInformativeText("Were corners rightly located?\n");
        msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox1.setDefaultButton(QMessageBox::Yes);
        int ret = msgBox1.exec();
        if(ret==QMessageBox::Yes){
            printf("OK! Found %d corners\n",numSquares);
            if(corners[0].y > corners[numSquares-1].y ){
                printf("found corners[0].y > corners[numSquares-1].y!!!\n");
                for(int iCor=0;iCor<numSquares;iCor++){
                    px[Np+numSquares-1-iCor][0]=corners[iCor].x;
                    px[Np+numSquares-1-iCor][1]=corners[iCor].y;
                }
            }
            else{
                for(int iCor=0;iCor<numSquares;iCor++){
                    px[Np+iCor][0]=corners[iCor].x;
                    px[Np+iCor][1]=corners[iCor].y;
                }
            }
        }
        else if(ret==QMessageBox::No)
            return;
    }
    imshow("Camera1",img1);
    imshow("Camera2",img2);
    waitKey(100);
    Np=Np*2;
    printf("Chessboard Ncorn=%d: Npoints = %d  Ntargets = %d\n",Ncorn,Np,2*Nt);
    fflush(stdout);
    if(2*Nt!=Np)
        printf("Warning: Npoints < Ntargets !!!!\n");

    Mat img3,img4,im3_with_keypoints,im4_with_keypoints;
    if(iSetCamMeth<2){
        PathPro="";
        PathPro=QFileDialog::getExistingDirectory(this, tr("Open Directory for setPointSource"),
                                                  dir+"/setSource",
                                                  QFileDialog::ShowDirsOnly);
        if(PathPro.isEmpty()){
            printf("abort!\n");
            return;
        }
        PathPro=PathPro+"/";
        zMirror=ui->dSB_zMirror->value();
        NbCam1=0;
        NbCam2=0;
        Mat imgCam3=imread((PathPro+"imgCam1.bmp").toStdString(),IMREAD_ANYDEPTH);//IMREAD_UNCHANGED);
        cout << (PathPro+"imgCam1.bmp").toStdString()<<"\n";
        if(!imgCam3.data){
            cout << (PathPro+"imgCam1.bmp").toStdString() <<" for viewImg is not accessible!!!\n"<<"\n";
            return;
        }
        Mat imgCam4=imread((PathPro+"imgCam2.bmp").toStdString(),IMREAD_ANYDEPTH);//IMREAD_UNCHANGED);
        cout << (PathPro+"imgCam2.bmp").toStdString() << "\n";
        if(!imgCam4.data){
            cout << (PathPro+"imgCam2.bmp").toStdString() <<" for viewImg is not accessible!!!\n"<<"\n";
            return;
        }
        state=ui-> checkBox_undist ->checkState();
        if(state==Qt::Checked)
            undistort(imgCam3,img3,cameraMatrix_1,distCoeffs_1);
        else
            img3=imgCam3;
        if(state==Qt::Checked)
            undistort(imgCam4,img4,cameraMatrix_2,distCoeffs_2);
        else
            img4=imgCam4;
        // set up the parameters (check the defaults in opencv's code in blobdetector.cpp)
        //cv::SimpleBlobDetector::Params params;

        double value;
        // Filter by treshold
        params.thresholdStep=5;
        params.minThreshold= ui -> sB_gL -> value();
        params.maxThreshold=255;

        //Filter by distance
        minDis=ui->doubleSpinBox_minDist->value();
        params.minDistBetweenBlobs = float(minDis);

        // Filter by Area
        state=ui->checkBox_Area->checkState();
        if(state==Qt::Checked){
            params.filterByArea = true;
            value=ui->doubleSpinBox_minArea->value();
            params.minArea = float(value);
            params.maxArea = 500.0f;
        }
        else
            params.filterByArea = false;

        // Filter by Circularity
        state=ui->checkBox_circularity->checkState();
        if(state==Qt::Checked){
            params.filterByCircularity = true;
            value=ui->doubleSpinBox_minCircularity->value();
            params.minCircularity =float(value);
        }
        else
            params.filterByCircularity = false;

        //filter by Convexity
        state=ui->checkBox_convexity->checkState();
        if(state==Qt::Checked){
            params.filterByConvexity = true;
            value=ui->doubleSpinBox_minConvexity->value();
            params.minConvexity =float(value);
        }
        else
            params.filterByConvexity = false;

        // Filter by Inertia
        state=ui->checkBox_inertia->checkState();
        if(state==Qt::Checked){
            params.filterByInertia = true;
            value=ui->doubleSpinBox_minInertia->value();
            params.minInertiaRatio =float(value);
        }
        else
            params.filterByInertia = false;

        params.filterByColor = false;
        // Filter by Circularity
        params.filterByCircularity = false;
        // ... any other params you don't want default value

        // set up and create the detector using the parameters
        //cv::SimpleBlobDetector blob_detector(params);
        cv::Ptr<cv::SimpleBlobDetector> blob_detector = cv::SimpleBlobDetector::create(params);

        for(int iCam=1;iCam<3;iCam++){
            // detect!
            vector<cv::KeyPoint> keypoints;

            // Draw detected blobs as red circles.
            // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures the size of the circle corresponds to the size of blob
            if(iCam==1){
                blob_detector -> detect(img3, keypoints);
                drawKeypoints(img3, keypoints, im3_with_keypoints, Scalar(0,255,0), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );
                imshow("ROIcam1", im3_with_keypoints );
            }
            else{
                blob_detector -> detect(img4, keypoints);
                drawKeypoints(img4, keypoints, im4_with_keypoints, Scalar(0,255,0), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );
                imshow("ROIcam2", im4_with_keypoints );
            }
            msgBox1.setText("ATTENTION!");
            msgBox1.setInformativeText("Were blobs rightly located?\n");
            msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
            msgBox1.setDefaultButton(QMessageBox::Yes);
            int ret = msgBox1.exec();
            if(ret==QMessageBox::No)
                return;
            // extract the x y coordinates of the keypoints:
            int nBlobs = keypoints.size();
            printf("Found nBlobs = %d\n",nBlobs);
            if(nBlobs>0){
                int iw=0;
                double blobRes[nBlobs][3];
                for (int i=0; i<nBlobs; i++){
                    float J = keypoints[i].pt.x;
                    float I = keypoints[i].pt.y;
                    float S = keypoints[i].size;
                    iw=i;
                    if(i>0){
                        while(J<blobRes[iw-1][0]){
                            iw--;
                        }
                        for(int ii=i-1;ii>=iw;ii--){
                            blobRes[ii+1][0]=blobRes[ii][0];
                            blobRes[ii+1][1]=blobRes[ii][1];
                            blobRes[ii+1][2]=blobRes[ii][2];
                        }
                    }
                    blobRes[iw][0]=J;
                    blobRes[iw][1]=I;
                    blobRes[iw][2]=S;
                }
                if(iCam==1)
                    NbCam1=nBlobs;
                else
                    NbCam2=nBlobs;
                int ibC100=-1,iB=0,iW=0;
                double Dj,Di,Djsum=0.,DjMax=-1000.,DjMin=1000.;
                //printf("iBlob=0 j=%f i=%f\n",blobRes[0][0],blobRes[0][1]);
                for (int ib=0; ib<nBlobs-1; ib++){
                    Dj=blobRes[iB+1][0]-blobRes[iB][0];
                    Di=blobRes[iB+1][1]-blobRes[iB][1];
                    if(sqrt(Dj*Dj+Di*Di)<minDis){//merge the two blobs
                        blobRes[iB][0]=0.5*(blobRes[iB][0]+blobRes[iB+1][0]);
                        blobRes[iB][1]=0.5*(blobRes[iB][1]+blobRes[iB+1][1]);
                        for(int k=iB+1;k<nBlobs-1;k++){//move 1 step back
                            for(int kk=0;kk<3;kk++)
                                blobRes[k][kk]=blobRes[k+1][kk];
                        }
                        nBlobs--;
                        //printf("two blobs were merged: %d with %d\nNew nBlobs=%d\n",iB,iB+1,nBlobs);
                        continue;
                    }
                    //printf("iBlob = %d j=%f i=%f S=%f Di= %f Dj= %f\n",
                    //       iB+1,blobRes[iB+1][0],blobRes[iB+1][1],blobRes[iB+1][2],Di,Dj);
                    if(ibC100==-1 && Di<-5. && (blobRes[iB+2][1]-blobRes[iB+1][1])>10){
                        ibC100=iB+1;//iCenter== point source N. 100
                        //printf("found ibC100=%d\n",ibC100);
                    }
                    else if(ibC100==-1 && Di>5. && (blobRes[iB+2][1]-blobRes[iB+1][1])<-10){
                        double c100[3];
                        for(int k=0;k<3;k++)
                            c100[k]=blobRes[iB+2][k];
                        for(int k=0;k<3;k++)
                            blobRes[iB+2][k]=blobRes[iB+1][k];
                        for(int k=0;k<3;k++)
                            blobRes[iB+1][k]=c100[k];
                        ibC100=iB+1;//iCenter== point source N. 100
                        //printf("found ibC100=%d\n",ibC100);
                        //printf("\tiBlob = %d j=%f i=%f\n",
                        //       iB+1,blobRes[iB+1][0],blobRes[iB+1][1]);
                        //printf("\tiBlob = %d j=%f i=%f\n",
                        //       iB+2,blobRes[iB+2][0],blobRes[iB+2][1]);
                    }
                    else if(abs(Di)<5.){
                        DjMin=min(DjMin,Dj);
                        DjMax=max(DjMax,Dj);
                        Djsum=Djsum+Dj;
                    }
                    iB++;
                }
                Dj=Djsum/double(nBlobs-2);
                //                    if(DjMax>1.75*Dj || DjMin<0.4*Dj)
                //                        iW=1;
                //                    else
                //                        iW=0;
                printf("<Dj>=%f\tiCenter= %d DjMin=%f DjMax=%f\n",
                       Dj,ibC100,DjMin,DjMax);
                if(ibC100==-1){
                    printf("abort: iCenter NOT found!!!\n");
                    fflush(stdout);
                    return;
                }
                else{
                    int nSource;
                    int ib;
                    int iStep=-1;
                    for (int i=0; i<nBlobs; i++){
                        ib=ibC100-i;
                        if(ib<0){
                            ib=ibC100+(i-ibC100);
                            iStep=1;
                        }
                        if(ib == ibC100 || ib == ibC100+1){
                            if(blobRes[ib][1] > blobRes[ibC100][1])
                                nSource=100;
                            else
                                nSource=99;
                        }
                        if(iW==1)
                            printf("i=%d ib=%d ",
                                   i,ib);
                        if(ib-iStep>=0 && ib-iStep<nBlobs && abs(blobRes[ib][0]-blobRes[ib-iStep][0]) > 1.75*abs(Dj)){
                            int DnS=NINT(abs((blobRes[ib][0]-blobRes[ib-iStep][0])/Dj));
                            nSource=nSource+DnS*double(-iStep);
                        }
                        else
                            nSource=nSource-iStep;
                        if(iCam==1){
                            Bcam1[ib][0]=blobRes[ib][0];
                            Bcam1[ib][1]=blobRes[ib][1];
                            Bcam1[ib][2]=nSource;
                            px[2*Nt+ib][0]=Bcam1[ib][0];
                            px[2*Nt+ib][1]=Bcam1[ib][1];
                            //printf("Bcam1[%d]:\tj= %f\ti= %f\tnSource= %d\n",
                            //       ib,Bcam1[ib][0],Bcam1[ib][1],NINT(Bcam1[ib][2]));
                        }
                        else{
                            Bcam2[ib][0]=blobRes[ib][0];
                            Bcam2[ib][1]=blobRes[ib][1];
                            Bcam2[ib][2]=nSource;
                            px[2*Nt+NbCam1+ib][0]=Bcam2[ib][0];
                            px[2*Nt+NbCam1+ib][1]=Bcam2[ib][1];
                            //printf("Bcam2[%d]:\tj= %f\ti= %f\tnSource= %d\n",
                            //       ib,Bcam2[ib][0],Bcam2[ib][1],NINT(Bcam2[ib][2]));
                        }
                    }
                }
            }
        }
        printf("Source reflection: NbCam1= %d\tNbCam2= %d\n",NbCam1,NbCam2);
        if(NbCam1==0 || NbCam2==0){
            printf("process abborted!\n");
            return;
        }
    }
    else{
        NbCam1=0;
        NbCam2=0;
    }
    msgBox.setText("ATTENTION!");
    msgBox.setInformativeText("Do you want to launch the best fit?\n");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::Yes);
    ret = msgBox.exec();
    if(ret==QMessageBox::No)
        return;
    // impostazione e lancio fit
    Pc_1[0]= ui-> dSB_Xc_1 -> value();
    Pc_1[1]= ui-> dSB_Yc_1 -> value();
    Pc_1[2]= ui-> dSB_Zc_1 -> value();
    EA_1[0]= ui-> dSB_Yaw_1 -> value(); //yaw
    EA_1[1]= ui-> dSB_Pitch_1 -> value(); //pitch
    EA_1[2]= ui-> dSB_Roll_1 -> value(); //roll
    printf("Initial values: Pc_1= %f\t%f\t%f\tEA_1= %f\t%f\t%f\n",Pc_1[0],Pc_1[1],Pc_1[2],EA_1[0],EA_1[1],EA_1[2]);
    for(int i=0;i<3;i++)
        EA_1[i]= EA_1[i]/180.*pig;
    Dx21=ui->dSB_Dx21->value();
    Dz21=ui->dSB_Dz21->value();
    //Pc_2[0]= ui-> dSB_Xc_2 -> value();
    Pc_2[0]=Pc_1[0]+Dx21;
    Pc_2[1]= ui-> dSB_Yc_2 -> value();
    //Pc_2[2]= ui-> dSB_Zc_2 -> value();
    Pc_2[2]=Pc_1[2]+Dz21;
    EA_2[0]= ui-> dSB_Yaw_2 -> value(); //yaw
    EA_2[1]= ui-> dSB_Pitch_2 -> value(); //pitch
    EA_2[2]= ui-> dSB_Roll_2 -> value(); //roll
    printf("Initial values: Pc_2= %f\t%f\t%f\tEA_2= %f\t%f\t%f\n",Pc_2[0],Pc_2[1],Pc_2[2],EA_2[0],EA_2[1],EA_2[2]);

    for(int i=0;i<3;i++)
        EA_2[i]= EA_2[i]/180.*pig;
    fx_1= ui-> lineEdit_fx_1 ->text().toDouble();
    fy_1= ui-> lineEdit_fy_1 ->text().toDouble();
    fx_2= ui-> lineEdit_fx_2 ->text().toDouble();
    fy_2= ui-> lineEdit_fy_2 ->text().toDouble();
    cx_1= ui-> lineEdit_cx_1 ->text().toDouble();
    cy_1= ui-> lineEdit_cy_1 ->text().toDouble();
    cx_2= ui-> lineEdit_cx_2 ->text().toDouble();
    cy_2= ui-> lineEdit_cy_2 ->text().toDouble();
    cameraMatrix_1.at<double>(0,0)=fx_1;
    cameraMatrix_1.at<double>(1,1)=fy_1;
    cameraMatrix_1.at<double>(0,2)=cx_1;
    cameraMatrix_1.at<double>(1,2)=cy_1;
    cameraMatrix_2.at<double>(0,0)=fx_2;
    cameraMatrix_2.at<double>(1,1)=fy_2;
    cameraMatrix_2.at<double>(0,2)=cx_2;
    cameraMatrix_2.at<double>(1,2)=cy_2;
    fcam_1=(fx_1*pxdimX+fy_1*pxdimY)/2.;
    fcam_2=(fx_2*pxdimX+fy_2*pxdimY)/2.;
    ui-> dSB_fcam_1  -> setValue(fcam_1);
    ui-> dSB_fcam_2  -> setValue(fcam_2);
    Ps[0]= ui-> dSB_Xs -> value();
    Ps[1]= ui-> dSB_Ys -> value();
    Ps[2]= ui-> dSB_Zs -> value();
    Es[0]= ui-> dSB_Yaw_s -> value(); //yaw
    Es[1]= ui-> dSB_Pitch_s -> value(); //pitch
    Es[2]= ui-> dSB_Roll_s -> value(); //roll
    printf("Initial values: Ps= %f\t%f\t%f\tEs= %f\t%f\t%f\n",Ps[0],Ps[1],Ps[2],Es[0],Es[1],Es[2]);
    for(int i=0;i<3;i++)
        Es[i]= Es[i]/180.*pig;

    int info=kernelSetCam(iSetCamMeth);

    cameraMatrix_1.at<double>(0,0)=fx_1;
    cameraMatrix_1.at<double>(1,1)=fy_1;
    cameraMatrix_1.at<double>(0,2)=cx_1;
    cameraMatrix_1.at<double>(1,2)=cy_1;
    cameraMatrix_2.at<double>(0,0)=fx_2;
    cameraMatrix_2.at<double>(1,1)=fy_2;
    cameraMatrix_2.at<double>(0,2)=cx_2;
    cameraMatrix_2.at<double>(1,2)=cy_2;
    fcam_1=(fx_1*pxdimX+fy_1*pxdimY)/2.;
    fcam_2=(fx_2*pxdimX+fy_2*pxdimY)/2.;

    printf("sqrt(chi2fin)= %f\n",sqrt(chi2fin));
    ui->lineEdit_chi2->setText(QString::number(chi2fin/Nt/2.));
    double DRMS=sqrt(chi2fin/Nt/2.);
    ui-> dSB_Xc_1 -> setValue(Pc_1[0]);
    ui-> dSB_Yc_1 -> setValue(Pc_1[1]);
    ui-> dSB_Zc_1 -> setValue(Pc_1[2]);
    ui-> dSB_Yaw_1   -> setValue(EA_1[0]*180./pig); //yaw
    ui-> dSB_Pitch_1 -> setValue(EA_1[1]*180./pig); //pitch
    ui-> dSB_Roll_1  -> setValue(EA_1[2]*180./pig); //roll
    ui-> dSB_Xc_2 -> setValue(Pc_2[0]);
    ui-> dSB_Yc_2 -> setValue(Pc_2[1]);
    ui-> dSB_Zc_2 -> setValue(Pc_2[2]);
    ui-> dSB_Yaw_2   -> setValue(EA_2[0]*180./pig); //yaw
    ui-> dSB_Pitch_2 -> setValue(EA_2[1]*180./pig); //pitch
    ui-> dSB_Roll_2  -> setValue(EA_2[2]*180./pig); //roll
    ui-> lineEdit_fx_1 ->setText(QString::number(fx_1));
    ui-> lineEdit_fy_1 ->setText(QString::number(fy_1));
    ui-> lineEdit_fx_2 ->setText(QString::number(fx_2));
    ui-> lineEdit_fy_2 ->setText(QString::number(fy_2));
    ui-> lineEdit_cx_1 ->setText(QString::number(cx_1));
    ui-> lineEdit_cy_1 ->setText(QString::number(cy_1));
    ui-> lineEdit_cx_2 ->setText(QString::number(cx_2));
    ui-> lineEdit_cy_2 ->setText(QString::number(cy_2));
//    ui->lineEdit_k1_1->setText(QString::number(distCoeffs_1.at<double>(0)));
//    ui->lineEdit_k2_1->setText(QString::number(distCoeffs_1.at<double>(1)));
//    ui->lineEdit_p1_1->setText(QString::number(distCoeffs_1.at<double>(2)));
//    ui->lineEdit_p2_1->setText(QString::number(distCoeffs_1.at<double>(3)));
//    ui->lineEdit_k3_1->setText(QString::number(distCoeffs_1.at<double>(4)));
//    ui->lineEdit_k1_2->setText(QString::number(distCoeffs_2.at<double>(0)));
//    ui->lineEdit_k2_2->setText(QString::number(distCoeffs_2.at<double>(1)));
//    ui->lineEdit_p1_2->setText(QString::number(distCoeffs_2.at<double>(2)));
//    ui->lineEdit_p2_2->setText(QString::number(distCoeffs_2.at<double>(3)));
//    ui->lineEdit_k3_2->setText(QString::number(distCoeffs_2.at<double>(4)));
    ui-> dSB_fcam_1  -> setValue(fcam_1);
    ui-> dSB_fcam_2  -> setValue(fcam_2);
    ui-> dSB_Xs    -> setValue(Ps[0]);
    ui-> dSB_Ys    -> setValue(Ps[1]);
    ui-> dSB_Zs    -> setValue(Ps[2]);
    ui-> dSB_Yaw_s   -> setValue(Es[0]*180./pig); //yaw
    ui-> dSB_Pitch_s -> setValue(Es[1]*180./pig); //pitch
    ui-> dSB_Roll_s  -> setValue(Es[2]*180./pig); //roll
    ui-> lineEdit_chi2 ->setText(QString::number(chi2fin));
    if(info>3)
        printf("Fit: info = %d\n\tchi2 = %f Drms = %f\n",info,chi2fin,DRMS);
    //traccia dot sui target e sul punto calcolato
    int iMin,iMax;
    Mat Cimg;
    x4set=ui->dSB_x4set->value();
    iXlab=x4set-(Pc_1[0]+Pc_2[0])/2.;
    double errXcur=abs(Pc_1[0]-Pc_2[0])/2.;
    ui->lineEdit_xLAB->setText(QString::number(iXlab)+" +- "+QString::number(errXcur));
    XlabCam1=Pc_1[0];
    ui->lineEdit_XlabCam1->setText(QString::number(XlabCam1));
    XlabCam2=Pc_2[0];
    ui->lineEdit_XlabCam2->setText(QString::number(XlabCam2));
    for(int iCam=1;iCam<3;iCam++){
        if(iCam==1){
            Cimg=img1;
            iMin=0;
            iMax=Nt;
        }
        else{
            Cimg=img2;
            iMin=Nt;
            iMax=2*Nt;
        }

        for(int i=iMin;i<iMax;i++){
            if (px[i][0] < 0.9e+8){
                Pt1.x=NINT(px[i][0]);
                Pt1.y=NINT(px[i][1]);
                //circle(Cimg, Pt1, 10, Scalar(65535), -1);
                circle(Cimg, Pt1, 10, Scalar(255), -1);
                Pt1.x=NINT(px[i][2]-10.);
                Pt1.y=NINT(px[i][3]-10);
                Pt2.x=Pt1.x+20;
                Pt2.y=Pt1.y+20;
                //rectangle(Cimg, Pt1,Pt2, Scalar(65535),1, 8, 0 );
                rectangle(Cimg, Pt1,Pt2,Scalar(255),1, 8, 0 );
            }
        }
        if(iCam==1)
            imshow("Camera1",Cimg);
        else
            imshow("Camera2",Cimg);
        QString saveimg=dir+"/img"+QString::number(iCam)+"WithTarget.JPEG";
        imwrite(saveimg.toStdString(),Cimg);
        waitKey(10);
    }
    if(iSetCamMeth<2){
    for(int i=2*Nt;i<2*Nt+NbCam1;i++){
        Pt1.x=NINT(px[i][2]-1.);
        Pt1.y=NINT(px[i][3]-1.);
        Pt2.x=Pt1.x+2;
        Pt2.y=Pt1.y+2;
        rectangle(im3_with_keypoints, Pt1,Pt2,Scalar(255),1, 8, 0 );
    }
    imshow("ROIcam1",im3_with_keypoints);
    imwrite(dir.toStdString()+"/image1WithBlob.JPEG",im3_with_keypoints);
    for(int i=2*Nt+NbCam1;i<2*Nt+NbCam1+NbCam2;i++){
        Pt1.x=NINT(px[i][2]-1.);
        Pt1.y=NINT(px[i][3]-1.);
        Pt2.x=Pt1.x+2;
        Pt2.y=Pt1.y+2;
        rectangle(im4_with_keypoints, Pt1,Pt2,Scalar(255),1, 8, 0 );
    }
    imshow("ROIcam2",im4_with_keypoints);
    imwrite(dir.toStdString()+"/image2WithBlob.JPEG",im4_with_keypoints);
    }
    waitKey(10);
    img1.release();
    img2.release();
    img3.release();
    img4.release();
    imgCam1.release();
    imgCam2.release();
    im3_with_keypoints.release();
    im4_with_keypoints.release();
    msgBox1.setText("ATTENTION!");
    msgBox1.setInformativeText("Do you want to save the new values?\n");
    msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox1.setDefaultButton(QMessageBox::Yes);
    int ret2 = msgBox1.exec();
    if(ret2==QMessageBox::Yes)
        callSaveAcam();
    else
        readAcam(fileAcam);

    QFile File(dir+"/PxPyTargets.txt");
    if (!File.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&File);
    out<<"Ncamera\tNpoint\tx_exp\ty_exp\tx_cal\ty_cal\n";
    for(int i=0;i<2*Nt+NbCam1+NbCam2;i++){
        if(i<Nt)
            out<<"cb1\t"<<i;
        else if(i<2*Nt)
            out<<"cb2\t"<<i-Nt;
        else if(i<2*Nt+NbCam1)
            out<<"cam1\t"<<i-2*Nt;
        else
            out<<"cam2\t"<<i-2*Nt-NbCam1;
        for(int ii=0;ii<4;ii++)
            out << "\t" << px[i][ii];
        out << "\n";
    }
    File.close();
}





void VISproPT::scan(){
    double newET=ui->dSBexpTime -> value();
    QString command;
    if(NINT(newET)!=NINT(expTime)){
        command=part1cmd+"set_expo_cam1.cmd "+QString::number(int(newET));
        system(command.toStdString().c_str());
        printf("set expTime_cam1 to %d\n",int(newET));
        command=part1cmd+"set_expo_cam2.cmd "+QString::number(int(newET));
        system(command.toStdString().c_str());
        printf("set expTime_cam2 to =%d\n",int(newET));
        expTime=newET;
    }
    ui->lineEdit_status ->setText("scan...");
    printf("-> launch run_cycle\n");
    Xmin=ui->dSB_Xmin->value();
    Xmax=ui->dSB_Xmax->value();
    Xstep=ui->dSB_Xstep->value();
    command=part1cmd+"run_cycle.cmd\t"+
            QString::number(Xmin)+"\t"+QString::number(Xmax)+"\t"+QString::number(Xstep);
    system(command.toStdString().c_str());
    ui->lineEdit_status ->setText("ready");
}


void VISproPT::process(){
//    QMessageBox msgBox1;
//    msgBox1.setText("ATTENTION!");
//    msgBox1.setInformativeText("Before processing you must have already load the right setting! It is the case?");
//    msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
//    msgBox1.setDefaultButton(QMessageBox::Yes);
//    int ret = msgBox1.exec();
//    if(ret==QMessageBox::No)
//        return;
    readParam();
    ui->lineEdit_status ->setText("process....");
    Qt::CheckState state;
    state=ui->checkBox_artCor->checkState();
    iArtCor=0;
    if(state==Qt::Checked)
        iArtCor=1;
    int nFrame;
    double value;
    iWr=ui->spinBox_iWr->value();
    //Mat img,imgc,im_with_keypoints;
    //namedWindow("image", WINDOW_AUTOSIZE );
    dirMeas=ui->lineEdit_img->text();
    if(dirMeas.isEmpty()){
        iTot=-1;//initialization
        dirMeas = QFileDialog::getExistingDirectory(this, tr("Open Directory for process"),
                                                    dir,
                                                    QFileDialog::ShowDirsOnly);
        printf("dirMeas= %s\n",dirMeas.toStdString().c_str());
        ui->lineEdit_img->setText(dirMeas);
        Xmax=ui->dSB_Xmax->value();
        Xmin=ui->dSB_Xmin->value();
        Xstep=ui->dSB_Xstep->value();
        Xstart=Xmax;
        nFrame=NINT((Xmax-Xmin)/Xstep)+1;
        Xstep=-Xstep;
        x4set=ui->dSB_x4set->value();
        printf("nFrame_expected=%d xStep=%f\n",nFrame,Xstep);

        // set up the parameters (check the defaults in opencv's code in blobdetector.cpp)
        //cv::SimpleBlobDetector::Params params;

        // Filter by treshold
        params.thresholdStep=5;
        params.minThreshold= ui -> sB_gL -> value();
        params.maxThreshold=255;

        //Filter by distance
        minDis=ui->doubleSpinBox_minDist->value();
        params.minDistBetweenBlobs = float(minDis);

        // Filter by Area
        state=ui->checkBox_Area->checkState();
        if(state==Qt::Checked){
            params.filterByArea = true;
            value=ui->doubleSpinBox_minArea->value();
            params.minArea = float(value);
            params.maxArea = 500.0f;
        }
        else
            params.filterByArea = false;

        // Filter by Circularity
        state=ui->checkBox_circularity->checkState();
        if(state==Qt::Checked){
            params.filterByCircularity = true;
            value=ui->doubleSpinBox_minCircularity->value();
            params.minCircularity =float(value);
        }
        else
            params.filterByCircularity = false;

        //filter by Convexity
        state=ui->checkBox_convexity->checkState();
        if(state==Qt::Checked){
            params.filterByConvexity = true;
            value=ui->doubleSpinBox_minConvexity->value();
            params.minConvexity =float(value);
        }
        else
            params.filterByConvexity = false;

        // Filter by Inertia
        state=ui->checkBox_inertia->checkState();
        if(state==Qt::Checked){
            params.filterByInertia = true;
            value=ui->doubleSpinBox_minInertia->value();
            params.minInertiaRatio =float(value);
        }
        else
            params.filterByInertia = false;

        params.filterByColor = false;
        // Filter by Circularity
        params.filterByCircularity = false;
        //params.minCircularity = 0.8f;

        // ... any other params you don't want default value

        for(int icam=1;icam<=2;icam++){
            //printf("Camera #%d:   x,    y,   z = %f\t%f\t%f\n",icam,PA[0],PA[1],PA[2]);
            //printf("            yaw,pitch,roll = %f\t%f\t%f\n",PA[3],PA[4],PA[5]);
            //printf("            fx(px)= %f\tfy(px)= %f\tcx= %f\tcy= %f\n",fx,fy,cx,cy);
            //printf("Xstart = %f xStep = %f nFrame=%d \n",Xstart,xStep,nFrame);
            kernelpro(icam);
        }
    }
    printf("%d Ndata obtained by image processing\n",iTot);

    double Dz=0.,DzOld;
    double DELTAz=20.;
    double DzMaxThresh=ui->dSB_thresLoop->value();//(mm)
    state=ui-> checkBox_undist ->checkState();
    if(state==Qt::Checked)
        undst1=1;
    else
        undst1=0;
    if(state==Qt::Checked)
        undst2=1;
    else
        undst2=0;

    state=ui->checkBox_OptLoop->checkState();
    idealHight(1);
    int Nite=0;
    do{
        DzOld=Dz;
        slopeComputing();
        if(iTot<1)
            return;
        inpainting();
        shapeComputing();
        Dz=abs(Djiz[2]);
        DELTAz=abs(Dz-DzOld);
        printf("Dz=%f DzOld=%f DELTAz=%f\n",Dz,DzOld,DELTAz);
        printf("******* Larger hight difference = %f at AP[%d][%d]=(%f,%f)\n\n",
               Djiz[2],NINT(Djiz[0]),NINT(Djiz[1]),(jiAP[NINT(Djiz[0])][NINT(Djiz[1])][1]-iScenter)*DS,
                (jiAP[NINT(Djiz[0])][NINT(Djiz[1])][0]-jScenter)*DS);
        Nite++;
    }while(DELTAz > DzMaxThresh && state==Qt::Checked && Nite<100);
    res2GUI();

    //create countormap
    mapMat();
    ui->lineEdit_status ->setText("ready");
}


void VISproPT::shapeComputing(){
    printf("shapeComputing with jAP=%d iAP=%d\n",jAP,iAP);
    int jStep,iStep,j1,i1,j2,i2,jc,ic,jOK,iOK;
    double iTmp=0.;
    double tmp[Nj][Ni]={{0.}};
    double tmp2[Nj][Ni]={{0.}};

    //2D integration around each attaching point
    //printf("done!\n\n\t2D integration around each attaching point...");
    for(int iap=0;iap<iAP;iap++){
        for(int jap=0;jap<jAP;jap++){
            if(jiAP[jap][iap][2]>0){
                ic=jiAP[jap][iap][1];
                jc=jiAP[jap][iap][0];
                S[ jiAP[jap][iap][0] ][ jiAP[jap][iap][1] ][7]=S[ jiAP[jap][iap][0] ][ jiAP[jap][iap][1] ][3];
                //printf("z at the attaching point [%d][%d]=%f\n",jiAP[jap][iap][0], jiAP[jap][iap][1],S[ jiAP[jap][iap][0] ][ jiAP[jap][iap][1] ][3]);
                for(int j=jiAP[jap][iap][0];jsMIN<j;j--){
                    S[j-1][ic][7]=S[j][ic][7]-(S[j][ic][6]+S[j-1][ic][6])*DS/2.;
                    //printf("right range: %d -> %d hight= %f Dy= %f ideal= %f\n",j,j-1,S[j-1][ic][7],(S[j][ic][6]+S[j-1][ic][6])/2.,(S[j][ic][9]+S[j-1][ic][9])/2.);
                }
                for(int j=jiAP[jap][iap][0];j<jsMAX;j++){
                    S[j+1][ic][7]=S[j][ic][7]+(S[j][ic][6]+S[j+1][ic][6])*DS/2.;
                    //printf("left range: %d -> %d hight= %f Dy= %f ideal= %f\n",j,j+1,S[j+1][ic][7],(S[j][ic][6]+S[j+1][ic][6])/2.,(S[j][ic][9]+S[j+1][ic][9])/2.);
                }
                for(int i=jiAP[jap][iap][1];isMIN<i;i--){
                    S[jc][i-1][7]=S[jc][i][7]-(S[jc][i][5]+S[jc][i-1][5])*DS/2.;
                    //printf("down range: %d -> %d hight= %f Dx= %f ideale= %f\n",i,i-1,S[jc][i-1][7],(S[jc][i][5]+S[jc][i-1][5])/2.,(S[jc][i][8]+S[jc][i-1][8])/2.);
                }
                for(int i=jiAP[jap][iap][1];i<isMAX;i++){
                    S[jc][i+1][7]=S[jc][i][7]+(S[jc][i][5]+S[jc][i+1][5])*DS/2.;
                    //printf("top range: %d -> %d hight= %f Dx= %f ideale= %f\n",i,i+1,S[jc][i+1][7],(S[jc][i][5]+S[jc][i+1][5])/2.,(S[jc][i][8]+S[jc][i+1][8])/2.);
                }
                for(int is=1;is<=4;is++){
                    jStep=pow(-1,is);
                    iStep=pow(-1,NINT(double(is)/2.));
                    j1=jiAP[jap][iap][0];
                    i1=jiAP[jap][iap][1];
                    if(jStep<0)
                        j2=jsMIN;
                    else
                        j2=jsMAX;
                    if(iStep<0)
                        i2=isMIN;
                    else
                        i2=isMAX;
                    ic=i1;
                    iOK=1;
                    //printf("\n\t\t(%d , %d) -> (%d , %d)\n",j1,i1,j2,i2);
                    while(iOK==1){
                        ic=ic+iStep;
                        jc=j1;
                        jOK=1;
                        while(jOK==1){
                            jc=jc+jStep;
                            if(jc>=jsMIN && jc<=jsMAX && ic>=isMIN && ic<=isMAX && S[jc][ic][4]>0.){
                                S[jc][ic][7]=(S[jc-jStep][ic][7]+(S[jc][ic][6]+S[jc-jStep][ic][6])*DS/2.*jStep)/2.;//3.;
                                S[jc][ic][7]=S[jc][ic][7]+(S[jc][ic-iStep][7]+(S[jc][ic][5]+S[jc][ic-iStep][5])*DS/2.*iStep)/2.;//3.;
                                //S[jc][ic][7]=S[jc][ic][7]+(S[jc-jStep][ic-iStep][7]+(S[jc][ic][5]+S[jc-jStep][ic-iStep][5])*iStep*DS/2.+(S[jc][ic][6]+S[jc-jStep][ic-iStep][6])*DS*jStep/2.)/3.;
                                //printf("*** z in (%d , %d)= %f from (%d ,%d) & (%d ,%d)\n",jc,ic,S[jc][ic][7],jc-jStep,ic,jc,ic-iStep);
                                fflush(stdout);
                            }
                            if(jc==j2 || jc<jsMIN || jc>jsMAX) jOK=0;
                        }
                        if(ic==i2 || ic<isMIN || ic>isMAX)iOK=0;
                    }
                }
                for(int ii=isMIN;ii<=isMAX;ii++){
                    for(int jj=jsMIN;jj<=jsMAX;jj++){
                        tmp[jj][ii]=tmp[jj][ii]+S[jj][ii][7];
                        tmp2[jj][ii]=tmp2[jj][ii]+S[jj][ii][7]*S[jj][ii][7];
                    }
                }
                iTmp++;
            }
        }
    }
    if(iTmp>0){
        for(int ii=isMIN;ii<=isMAX;ii++){
            for(int jj=jsMIN;jj<=jsMAX;jj++){
                S[jj][ii][7]=tmp[jj][ii]/iTmp;
                S[jj][ii][11]=sqrt(tmp2[jj][ii]/iTmp-pow(tmp[jj][ii]/iTmp,2.))/sqrt(iTmp);
            }
        }
        Djiz[2]=0.;
        double Dz;
        for(int iap=0;iap<iAP;iap++){
            for(int jap=0;jap<jAP;jap++){
                if(jiAP[jap][iap][2]>0){
                   Dz=abs(S[ jiAP[jap][iap][0] ][ jiAP[jap][iap][1] ][7]-
                          S[ jiAP[jap][iap][0] ][ jiAP[jap][iap][1] ][3]);
                   //printf("jAP=%d iAP%d Dz=%f\n",jiAP[jap][iap][0],jiAP[jap][iap][1],Dz);
                   if(Dz>Djiz[2]){
                       Djiz[0]=jap;
                       Djiz[1]=iap;
                       Djiz[2]=Dz;
                   }
                }
            }
        }
    }
    printf("\tshapeComputig done with iTmp=%d\n",NINT(iTmp));
    fflush(stdout);
    ui->lineEdit_status ->setText("ready");
}


void VISproPT::idealShape(){
    idealHight(1);
    res2GUI();
}


void VISproPT::shape(){
    shapeComputing();
    res2GUI();
}


void VISproPT::slope(){
    readParam();
    slopeComputing();
    inpainting();
    res2GUI();
}



void VISproPT::inpainting(){
    double MaxPopRow=0.;
    double MaxPopCol=0.;
    double soglia=ui->dSB_inpaThres->value();
    //computing row population
    for(int i=0; i<Ni; i++){
        S[0][i][10]=0.;
        for(int j=0; j<Nj; j++){
            if(NINT(S[j][i][4])>0)
                S[0][i][10]++;
        }
        MaxPopRow=max(MaxPopRow,S[0][i][10]);
    }
    //compuitng column population
    for(int j=0; j<Nj; j++){
        S[j][0][10]=0.;
        for(int i=0; i<Ni; i++){
            if(NINT(S[j][i][4])>0)
                S[j][0][10]++;
        }
        MaxPopCol=max(MaxPopCol,S[j][0][10]);
    }
    printf("inpainting...\n\tMaxPopRow=%d\tMaxPopCol=%d\n",int(MaxPopRow),int(MaxPopCol));
    isMIN=0;
    while(S[0][isMIN][10]/MaxPopRow<soglia)
        isMIN++;
    //printf("isMIN=%d with population=%d\n",isMIN,int(S[0][isMIN][10]));
    isMAX=Ni;
    while(S[0][isMAX][10]/MaxPopRow<soglia)
        isMAX--;
    //printf("isMAX=%d with population=%d\n",isMAX,int(S[0][isMAX][10]));
    jsMIN=0;
    while(S[jsMIN][0][10]/MaxPopCol<soglia)
        jsMIN++;
    //printf("jsMIN=%d with population=%d\n",jsMIN,int(S[jsMIN][0][10]));
    jsMAX=Nj;
    while(S[jsMAX][0][10]/MaxPopCol<soglia)
        jsMAX--;
    //printf("jsMAX=%d with population=%d\n",jsMAX,int(S[jsMAX][0][10]));

    int maxData=0;
    double sum=0.;
    double nCell=0.;
    Mat sampled=Mat::zeros(isMAX-isMIN+1,jsMAX-jsMIN+1,CV_8UC1);
    for(int i=isMIN; i<=isMAX; i++){
        for(int j=jsMIN; j<=jsMAX; j++){
            if(NINT(S[j][i][4])>0){
                sampled.at<uchar>(isMAX-i,jsMAX-j)=255;
                maxData=max(maxData,NINT(S[j][i][4]));
                sum=sum+S[j][i][4];
                nCell++;
                double mod=0.;
                for(int k=0;k<3;k++)//averaging
                    mod=mod+S[j][i][k]*S[j][i][k];
                for(int k=0;k<3;k++)
                    S[j][i][k]=S[j][i][k]/sqrt(mod);
                S[j][i][5]=-S[j][i][0]/S[j][i][2];//dz/dx (partial x-derivative)
                S[j][i][6]=-S[j][i][1]/S[j][i][2];//dz/dy (partial y-derivative)
                //printf("(%d, %d): z= %f Dx= %f Dy= %f\n",j,i,S[j][i][3],S[j][i][5],S[j][i][6]);
            }
        }
    }
    printf("S matrix is sampled in the region\t[jMax %d][iMax %d]\n"
           "\t\t\t\t\t[jMin %d][iMin %d]\n",jsMAX,isMAX,jsMIN,isMIN);
    //printf("jScenter=%d iScenter=%d\n",jScenter,iScenter);
    printf("N. max data in the single cell = %d\tmean_value = %f\n",maxData,sum/nCell);
    printf("Starting inpainting ...");
    namedWindow("sampledSurface", WINDOW_NORMAL);
    imshow("sampledSurface",sampled);
    imwrite((dir+"/sampledSurface.jpg").toStdString(),sampled);
    for(int i=isMIN; i<=isMAX; i++){
        for(int j=jsMIN; j<=jsMAX; j++){
            if(NINT(S[j][i][4])==0){
                if(j-1>=jsMIN && j+1<=jsMAX){
                    if(S[j-1][i][4]>0 && S[j+1][i][4]>0){
                        for(int k=0; k<3 ; k++)
                            S[j][i][k]=S[j][i][k]+(S[j-1][i][k]+S[j+1][i][k])/2.;
                        S[j][i][4]++;
                    }
                }
                if(j-1>=jsMIN && j+1<=jsMAX && i-1 >=isMIN && i+1<=isMAX){
                    if(S[j-1][i-1][4]>0 && S[j+1][i+1][4]>0){
                        for(int k=0; k<3 ; k++){
                            S[j][i][k]=S[j][i][k]+(S[j-1][i-1][k]+S[j+1][i+1][k])/2.;
                            S[j][i][k]=S[j][i][k]+(S[j-1][i+1][k]+S[j+1][i-1][k])/2.;
                        }
                        S[j][i][4]++;
                    }
                }
                if(i-1>=isMIN && i+1<=isMAX){
                    if(S[j][i-1][4]>0 && S[j][i+1][4]>0){
                        for(int k=0; k<3 ; k++)
                            S[j][i][k]=S[j][i][k]+(S[j][i-1][k]+S[j][i+1][k])/2.;
                        S[j][i][4]++;
                    }
                }
                if(NINT(S[j][i][4])>0){
                    double mod=0.;
                    for(int k=0;k<3;k++)//averaging
                        mod=mod+S[j][i][k]*S[j][i][k];
                    for(int k=0; k<3 ; k++)
                        S[j][i][k]=S[j][i][k]/sqrt(mod);
                    S[j][i][5]=-S[j][i][0]/S[j][i][2];//dz/dx (partial x-derivative)
                    S[j][i][6]=-S[j][i][1]/S[j][i][2];//dz/dy (partial y-derivative)
                    //printf("cell(%d,%d) was inpainted by %d independent evaluations z=%f Dx= %f Dy= %f\n",j,i,NINT(S[j][i][4]),S[j][i][3],S[j][i][5],S[j][i][6]);
                }
                if(NINT(S[j][i][4])>0){
                    sampled.at<uchar>(isMAX-i,jsMAX-j)=255;
                    double mod=0.;
                    for(int k=0; k<3 ; k++)
                        mod=mod+S[j][i][k]*S[j][i][k];
                    if(mod>1.00001)
                        printf("warning mod=%f @ j=%d i=%d\n",mod,j,i);

                }
            }
        }
    }
    namedWindow("inpaintedSurface", WINDOW_NORMAL);
    imshow("inpaintedSurface",sampled);
    printf("\tinpainting done!\n");
    fflush(stdout);
}


double VISproPT::slopeChi2(){
    idealHight(0);
    nScellFill=0;
    double chi2=0.;
    for(int i=isMIN;i<=isMAX;i++){
        for(int j=jsMIN;j<=jsMAX;j++){
            if(NINT(S[j][i][4])>0){
                nScellFill++;
                chi2=chi2+(S[j][i][5]*1000.)*(S[j][i][5]*1000.);
                chi2=chi2+(S[j][i][6]*1000.)*(S[j][i][6]*1000.);
            }
        }
    }
    chi2=chi2/double(nScellFill)/2.;
    return(chi2);
}



void VISproPT::readParam(){
    printf("readParam\n");
    ui->lineEdit_status ->setText("busy");
    waitKey(10);
    Pc_1[0]= ui-> dSB_Xc_1 -> value();
    Pc_1[1]= ui-> dSB_Yc_1 -> value();
    Pc_1[2]= ui-> dSB_Zc_1 -> value();
    EA_1[0]= ui-> dSB_Yaw_1 -> value(); //yaw
    EA_1[1]= ui-> dSB_Pitch_1 -> value(); //pitch
    EA_1[2]= ui-> dSB_Roll_1 -> value(); //roll
    for(int i=0;i<3;i++)
        EA_1[i]= EA_1[i]/180.*pig;
    Pc_2[0]= ui-> dSB_Xc_2 -> value();
    Pc_2[1]= ui-> dSB_Yc_2 -> value();
    Pc_2[2]= ui-> dSB_Zc_2 -> value();
    EA_2[0]= ui-> dSB_Yaw_2 -> value(); //yaw
    EA_2[1]= ui-> dSB_Pitch_2 -> value(); //pitch
    EA_2[2]= ui-> dSB_Roll_2 -> value(); //roll
    for(int i=0;i<3;i++)
        EA_2[i]= EA_2[i]/180.*pig;
    fx_1=ui->lineEdit_fx_1-> text().toDouble();
    fy_1=ui->lineEdit_fy_1-> text().toDouble();
    cx_1=ui->lineEdit_cx_1-> text().toDouble();
    cy_1=ui->lineEdit_cy_1-> text().toDouble();
    cameraMatrix_1.at<double>(0,0)=fx_1;
    cameraMatrix_1.at<double>(1,1)=fy_1;
    cameraMatrix_1.at<double>(0,2)=cx_1;
    cameraMatrix_1.at<double>(1,2)=cy_1;
    fcam_1=(fx_1*pxdimX+fy_1*pxdimY)/2.;
    distCoeffs_1.at<double>(0)= ui->lineEdit_k1_1->text().toDouble();
    distCoeffs_1.at<double>(1)= ui->lineEdit_k2_1->text().toDouble();
    distCoeffs_1.at<double>(2)= ui->lineEdit_p1_1->text().toDouble();
    distCoeffs_1.at<double>(3)= ui->lineEdit_p2_1->text().toDouble();
    distCoeffs_1.at<double>(4)= ui->lineEdit_k3_1->text().toDouble();
    distCoeffs_1.at<double>(5)= ui->lineEdit_k4_1->text().toDouble();
    distCoeffs_1.at<double>(6)= ui->lineEdit_k5_1->text().toDouble();
    distCoeffs_1.at<double>(7)= ui->lineEdit_k6_1->text().toDouble();
    distCoeffs_1.at<double>(8)= ui->lineEdit_s1_1->text().toDouble();
    distCoeffs_1.at<double>(9)= ui->lineEdit_s2_1->text().toDouble();
    distCoeffs_1.at<double>(10)=ui->lineEdit_s3_1->text().toDouble();
    distCoeffs_1.at<double>(11)=ui->lineEdit_s4_1->text().toDouble();
    distCoeffs_1.at<double>(12)=ui->lineEdit_tx_1->text().toDouble();
    distCoeffs_1.at<double>(13)=ui->lineEdit_ty_1->text().toDouble();
    ui-> dSB_fcam_1  -> setValue(fcam_1);
    fx_2=ui->lineEdit_fx_2-> text().toDouble();
    fy_2=ui->lineEdit_fy_2-> text().toDouble();
    cx_2=ui->lineEdit_cx_2-> text().toDouble();
    cy_2=ui->lineEdit_cy_2-> text().toDouble();
    cameraMatrix_2.at<double>(0,0)=fx_2;
    cameraMatrix_2.at<double>(1,1)=fy_2;
    cameraMatrix_2.at<double>(0,2)=cx_2;
    cameraMatrix_2.at<double>(1,2)=cy_2;
    fcam_2=(fx_2*pxdimX+fy_2*pxdimY)/2.;
    distCoeffs_2.at<double>(0)= ui->lineEdit_k1_2->text().toDouble();
    distCoeffs_2.at<double>(1)= ui->lineEdit_k2_2->text().toDouble();
    distCoeffs_2.at<double>(2)= ui->lineEdit_p1_2->text().toDouble();
    distCoeffs_2.at<double>(3)= ui->lineEdit_p2_2->text().toDouble();
    distCoeffs_2.at<double>(4)= ui->lineEdit_k3_2->text().toDouble();
    distCoeffs_2.at<double>(5)= ui->lineEdit_k4_2->text().toDouble();
    distCoeffs_2.at<double>(6)= ui->lineEdit_k5_2->text().toDouble();
    distCoeffs_2.at<double>(7)= ui->lineEdit_k6_2->text().toDouble();
    distCoeffs_2.at<double>(8)= ui->lineEdit_s1_2->text().toDouble();
    distCoeffs_2.at<double>(9)= ui->lineEdit_s2_2->text().toDouble();
    distCoeffs_2.at<double>(10)=ui->lineEdit_s3_2->text().toDouble();
    distCoeffs_2.at<double>(11)=ui->lineEdit_s4_2->text().toDouble();
    distCoeffs_2.at<double>(12)=ui->lineEdit_tx_2->text().toDouble();
    distCoeffs_2.at<double>(13)=ui->lineEdit_ty_2->text().toDouble();
    ui-> dSB_fcam_2  -> setValue(fcam_2);
    if(nScellFill!=0)
        res2GUI();
    ui->lineEdit_status ->setText("ready");
}


void VISproPT::rmTilt(){
    printf("remove tilt...\n");
    double corrDx=0.,corrDy=0.;
    MinMax[0][0]=1.e6;
    MinMax[0][1]=-1.e6;
    MinMax[1][0]=1.e6;
    MinMax[1][1]=-1.e6;
    MinMax[2][0]=1.e6;
    MinMax[2][1]=-1.e6;
    for(int i=isMIN;i<=isMAX;i++){
        for(int j=jsMIN;j<=jsMAX;j++){
            if(NINT(S[j][i][4])>0){
                corrDx=(((Dz[1][1]-Dz[0][1])-(Dz[1][0]-Dz[0][0]))/
                        double(jiAP[0][1][1]-jiAP[0][0][1])*(i-jiAP[0][0][1])+(Dz[1][0]-Dz[0][0]))/
                        (double(jiAP[1][0][0]-jiAP[0][0][0])*DS);
                corrDy=(((Dz[1][0]-Dz[1][1])-(Dz[0][0]-Dz[0][1]))/
                        double(jiAP[1][0][0]-jiAP[0][0][0])*(j-jiAP[0][0][0])+(Dz[0][0]-Dz[0][1]))/
                        (double(jiAP[0][1][1]-jiAP[0][0][1])*DS);
                S[j][i][5]=S[j][i][5]-corrDx;
                S[j][i][6]=S[j][i][6]-corrDy;
                MinMax[0][0]=min(MinMax[0][0],S[j][i][5]);//slopeX
                MinMax[1][0]=min(MinMax[1][0],S[j][i][6]);//slopeY
                MinMax[2][0]=min(MinMax[2][0],S[j][i][7]);
                MinMax[0][1]=max(MinMax[0][1],S[j][i][5]);
                MinMax[1][1]=max(MinMax[1][1],S[j][i][6]);
                MinMax[2][1]=max(MinMax[2][1],S[j][i][7]);
            }
        }
    }
    double RangeY=MinMax[0][1]-MinMax[0][0];//slopeY
    double RangeX=MinMax[1][1]-MinMax[1][0];//slopeX
    double RangeZ=MinMax[2][1]-MinMax[2][0];//z
    printf("rmTilt Range:\n\tDx = [%f , %f]\tRange=%f\n\tDy = [%f , %f]\tRange=%f\n\t z = [%f , %f]\tRange=%f\n",
           MinMax[0][0],MinMax[0][1],RangeX,MinMax[1][0],MinMax[1][1],RangeY,MinMax[2][0],MinMax[2][1],RangeZ);
    shapeComputing();
    res2GUI();
}

void VISproPT::mapMat(){
    cout<<"mapMat: creation and plot of slope Maps"<<"\n";
    doIntFat();

    // creation and plot of slope Maps
    Mat M(2,2, CV_8UC3, Scalar(0,0,255));
    Mat slopeXdevMap(isMAX-isMIN+1,jsMAX-jsMIN+1,CV_8UC3,Scalar(150,150,150));
    Mat slopeYdevMap(isMAX-isMIN+1,jsMAX-jsMIN+1,CV_8UC3,Scalar(150,150,150));
    Mat ZdevMap(isMAX-isMIN+1,jsMAX-jsMIN+1,CV_8UC3,Scalar(150,150,150));
    Mat intFatMap(isMAX-isMIN+1,jsMAX-jsMIN+1,CV_8UC3,Scalar(150,150,150));
    double val;
    double rangeDevSlope=ui->dSB_rangeSlopeDev->value();
    rangeDevSlope=rangeDevSlope/1000.*2;
    double rangeDevZ=ui->dSB_rangeZdev->value();
    rangeDevZ=rangeDevZ*2;
    double RangeDx=MinMax[0][1]-MinMax[0][0];//devSlope dz/dx
    double RangeDy=MinMax[1][1]-MinMax[1][0];//devSlope dz/dy
    double RangeZ=MinMax[2][1]-MinMax[2][0];//devZ
    double sx=0.,sxx=0.,sy=0.,syy=0.,sz=0.,szz=0.;
    int nDat=0;
    for(int i=isMIN;i<=isMAX;i++){
        for(int j=jsMIN;j<=jsMAX;j++){
            fflush(stdout);
            if(S[j][i][4]>0){
                nDat++;

                //slopeX
                val=(S[j][i][5]-S[j][i][8])/rangeDevSlope+0.5;
                pxColor(val);
                for(int k=0;k<3;k++)
                    slopeXdevMap.at<Vec3b>(isMAX-i,jsMAX-j)[k] = static_cast<uint8_t>(Vservice[k]);
                sx=sx+S[j][i][5]-S[j][i][8];
                sxx=sxx+pow(S[j][i][5]-S[j][i][8],2.);

                //slopeY
                val=(S[j][i][6]-S[j][i][9])/rangeDevSlope+0.5;
                pxColor(val);
                for(int k=0;k<3;k++)
                    slopeYdevMap.at<Vec3b>(isMAX-i,jsMAX-j)[k] = static_cast<uint8_t>(Vservice[k]);
                sy=sy+S[j][i][6]-S[j][i][9];
                syy=syy+pow(S[j][i][6]-S[j][i][9],2.);

                //Z
                val=(S[j][i][7]-S[j][i][3])/rangeDevZ+0.5;
                pxColor(val);
                for(int k=0;k<3;k++)
                    ZdevMap.at<Vec3b>(isMAX-i,jsMAX-j)[k] = static_cast<uint8_t>(Vservice[k]);
                sz=sz+S[j][i][7]-S[j][i][3];
                szz=szz+pow(S[j][i][7]-S[j][i][3],2.);

                //intFat
                if(S[j][i][10]>=0.9999){
                    for(int k=0;k<3;k++)
                        intFatMap.at<Vec3b>(isMAX-i,jsMAX-j)[k]=255;
                }
                else {
                    pxColor(S[j][i][10]);
                    for(int k=0;k<3;k++)
                        intFatMap.at<Vec3b>(isMAX-i,jsMAX-j)[k] = static_cast<uint8_t>(Vservice[k]);
                }

            }
        }
    }

    namedWindow("devSlopeX Map", WINDOW_NORMAL);
    imshow("devSlopeX Map",slopeXdevMap);
    imwrite(dir.toStdString()+"/expo/devSlopeXmap.JPEG",slopeXdevMap);
    namedWindow("devSlopeY Map", WINDOW_NORMAL);
    imshow("devSlopeY Map",slopeYdevMap);
    imwrite(dir.toStdString()+"/expo/devSlopeYmap.JPEG",slopeYdevMap);
    namedWindow("devZ Map", WINDOW_NORMAL);
    imshow("devZ Map",ZdevMap);
    imwrite(dir.toStdString()+"/expo/devZmap.JPEG",ZdevMap);
    namedWindow("intFat Map",WINDOW_NORMAL);
    imshow("intFat Map",intFatMap);
    imwrite(dir.toStdString()+"/expo/intFat.JPEG",intFatMap);

    printf("Range deviation (nDat=%d in S matrix obtained by iTot=%d row data):\n"
           "\tslopeX = [%f , %f]\tRange=%f\tMean=%f\trms=%e\n"
           "\tslopeY = [%f , %f]\tRange=%f\tMean=%f\trms=%e\n"
           "\tz      = [%f , %f]\tRange=%f\tMean=%f\trms=%e\n",nDat,iTot,
           MinMax[0][0],MinMax[0][1],RangeDx,sx/nDat,sqrt(sxx/nDat-pow(sx/nDat,2.)),
           MinMax[1][0],MinMax[1][1],RangeDy,sy/nDat,sqrt(syy/nDat-pow(sy/nDat,2.)),
           MinMax[2][0],MinMax[2][1],RangeZ ,sz/nDat,sqrt(szz/nDat-pow(sz/nDat,2.)));
    printf("x_range: [%f , %f]\ty_range [%f , %f]\n",(isMIN-iScenter)*DS,(isMAX-iScenter)*DS,(jsMIN-jScenter)*DS,(jsMAX-jScenter)*DS);
    slopeXdevMap.release();
    slopeYdevMap.release();
    ZdevMap.release();
    intFatMap.release();
    fflush(stdout);

    //Save results on file

    //MatDz matrix
    expo(fileMatDz,7);

    //MatZrms
    expo(fileMatRMSz,11);

    //MatXslope matrix
    expo(fileMatXslope,5);

    //MatYslope matrix
    expo(fileMatYslope,6);

    //intFat matrix
    expo(fileMatIntFat,10);
}


void VISproPT::doIntFat(){
    intFat();
    stat(10);
    ui -> lineEdit_intFat -> setText(QString::number(Mean));
    printf("<intFat> = %f\n",Mean);
}


void VISproPT::expo(QString file2expo,int iq){
    int iq0;
    double secT;
    Qt::CheckState state;
    state=ui-> checkBox_exportDev ->checkState();
    if(state==Qt::Checked){
        if(iq==5)//Dx
            iq0=8;
        else if(iq==6)//Dy
            iq0=9;
        else if(iq==7)//z
            iq0=3;
        else if(iq==10)//intFat
            iq0=-1;
        else if(iq==11)//Zrms
            iq0=-1;
        else
            return;
    }
    else
        iq0=-1;
    QFile FileExpo(file2expo);
    if (!FileExpo.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&FileExpo);
    for(int i=isMIN; i<=isMAX; i++){
        for(int j=jsMAX; j>=jsMIN; j--){//QtiPLot pone la cella [0][0] nell'origine del riferimento
            if(NINT(S[j][i][4])>0){
                if(iq0>0)
                    secT=S[j][i][iq0];
                else
                    secT=0.;
                out << S[j][i][iq]-secT;
            }
            else
                out << "NAN";
            if(j>jsMIN)
                out << "\t";
            else
                out << "\n";
        }
    }
    FileExpo.close();
}


void VISproPT::stat(int iq){
    int iq0;
    double sum=.0,sum2=0.,ndat=0.,secT;
    if(iq==5)//Dx
        iq0=8;
    else if(iq==6)//Dy
        iq0=9;
    else if(iq==7)//z
        iq0=3;
    else if(iq==10)//intFat
        iq0=-1;
    else
        return;
    for(int i=isMIN;i<=isMAX;i++){
        for(int j=jsMIN;j<=jsMAX;j++){
            if(NINT(S[j][i][4])>0){
                if(iq0>0)
                    secT=S[j][i][iq0];
                else
                    secT=0.;
                sum =sum +S[j][i][iq]-secT;
                sum2=sum2+(S[j][i][iq]-secT)*(S[j][i][iq]-secT);
                ndat++;
            }
        }
    }
    Mean=sum/ndat;
    standDev=sum2/ndat-Mean*Mean;
    if(standDev>0.)
        standDev=sqrt(standDev);
}



void VISproPT::res2GUI(){
    printf("resGUI...\n");
    // computing mean value in S[j][i][k] and focal_length in central belt
    double sx=0.,sxx=0.,sxDy=0.,sDy=0.,sy=0.,syy=.0,syDx=.0,sDx=0.;
    double isMiddle=double(isMAX+isMIN)/2.;
    //double isRange =double(isMAX-isMIN);
    double jsMiddle=double(jsMAX+jsMIN)/2.;
    //double jsRange =double(jsMAX-jsMIN);
    double Ndat=0.;
    nScellFill=0;
    double chi2=0.;
    MinMax[0][0]=1.e6;
    MinMax[0][1]=-1.e6;
    MinMax[1][0]=1.e6;
    MinMax[1][1]=-1.e6;
    MinMax[2][0]=1.e6;
    MinMax[2][1]=-1.e6;
    for(int i=isMIN;i<=isMAX;i++){
        for(int j=jsMIN;j<=jsMAX;j++){
            if(NINT(S[j][i][4])>0){
                nScellFill++;
                //                double sum=0.;
                //                for(int k=0;k<3;k++){//
                //                    sum=sum+S[j][i][k]*S[j][i][k];
                //                }
                //                sum=sqrt(sum);
                //                for(int k=0;k<3;k++)//normalization
                //                    S[j][i][k]=S[j][i][k]/sum;
                //                S[j][i][5]=S[j][i][1]/S[j][i][2];//dz/dy (partial y-derivative)
                //                S[j][i][6]=S[j][i][0]/S[j][i][2];//dz/dx (partial x-derivative)
                chi2=chi2+(S[j][i][5]*1000.)*(S[j][i][5]*1000.);
                chi2=chi2+(S[j][i][6]*1000.)*(S[j][i][6]*1000.);
                MinMax[0][0]=min(MinMax[0][0],S[j][i][5]-S[j][i][8]);
                MinMax[1][0]=min(MinMax[1][0],S[j][i][6]-S[j][i][9]);
                MinMax[2][0]=min(MinMax[2][0],S[j][i][7]-S[j][i][3]);
                MinMax[0][1]=max(MinMax[0][1],S[j][i][5]-S[j][i][8]);
                MinMax[1][1]=max(MinMax[1][1],S[j][i][6]-S[j][i][9]);
                MinMax[2][1]=max(MinMax[2][1],S[j][i][7]-S[j][i][3]);
                //if(abs((double(i)-isMiddle)/isRange)<0.25){
                    Ndat++;
                    sx=sx+(isMiddle-double(i))*DS;
                    sxx=sxx+(isMiddle-double(i))*DS*(isMiddle-double(i))*DS;
                    sxDy=sxDy+(isMiddle-double(i))*DS*S[j][i][6];
                    sDy=sDy+S[j][i][6];
                    sy=sy-(jsMiddle-double(j))*DS;
                    syy=syy+(jsMiddle-double(j))*DS*(jsMiddle-double(j))*DS;
                    syDx=syDx-(jsMiddle-double(j))*DS*S[j][i][5];
                    sDx=sDx+S[j][i][5];
                //}
            }
        }
    }
    chi2=chi2/double(nScellFill)/2.;
    for(int i=0;i<2;i++){
        for(int j=0; j<2; j++)
            if(jiAP[j][i][2]>0)
                Dz[j][i]=S[ jiAP[j][i][0] ][ jiAP[j][i][1] ][7] - S[ jiAP[j][i][0] ][ jiAP[j][i][1] ][3];
    }
    if(jiAP[0][0][2]>0)
        ui -> lineEdit_z0_0 -> setText(QString::number(Dz[0][0]));
    else
        ui -> lineEdit_z0_0 -> setText("NA");

    if(jiAP[1][0][2]>0)
        ui -> lineEdit_z1_0 -> setText(QString::number(Dz[1][0]));
    else
        ui -> lineEdit_z1_0 -> setText("NA");
    if(jiAP[0][1][2]>0)
        ui -> lineEdit_z0_1 -> setText(QString::number(Dz[0][1]));
    else
        ui -> lineEdit_z0_1 -> setText("NA");
    if(jiAP[1][1][2]>0)
        ui -> lineEdit_z1_1 -> setText(QString::number(Dz[1][1]));
    else
        ui -> lineEdit_z1_1 ->setText("NA");

    printf("z(0,0) z(1,0)= %f\t%f\n",Dz[0][0],Dz[1][0]);
    printf("z(0,1) z(1,1)= %f\t%f\n",Dz[0][1],Dz[1][1]);
    fflush(stdout);
}



void VISproPT::idealHight(int iWarning){
    //given the considered panel, here attaching points are set and ideal shape computed
    iWr=ui->spinBox_iWr->value();
    double xsiA1=0.,xsiA2=0.,etaA1=0.,etaA2=0.,Dx,Dy,dz0holder=0.;
    int iPanel= ui-> comboBox_panel -> currentIndex();
    if(iWarning==1)
        printf("compunting ideal hight panel type %d......\n"
               "Ni=%d iScenter=%d DS=%f\n",iPanel,Ni,iScenter,DS);

    if(iPanel==0){//ENEA inner panel
        focal=1810.;// focal lenght
        xsiA1=375.4;  //attaching points
        xsiA2=1307.9;
        etaA1=-11.5;
        etaA2=203.8;
        Dx=sqrt(pow(xsiA2-xsiA1,2.)+pow(etaA2-etaA1,2.))/2.;
        Dy=720./2.;
        //dz0holder=95.20;
    }
    else if(iPanel==1){//ENEA outer panel
        focal=1810.;// focal lenght
        xsiA1=1910.8;  //attaching points
        xsiA2=2718.4;
        etaA1=470.8;
        etaA2=984.3;
        Dx=sqrt(pow(xsiA2-xsiA1,2.)+pow(etaA2-etaA1,2.))/2.;
        Dy=720./2.;
        //dz0holder=95.20;
    }
    else if(iPanel==2){//RR inner panel
        focal=1710.;// focal lenght
        xsiA1=372.25;  //attaching points
        xsiA2=1325.52;
        etaA1=-0.76;
        etaA2=234.46;
        Dx=sqrt(pow(xsiA2-xsiA1,2.)+pow(etaA2-etaA1,2.))/2.;
        Dy=996./2.;
        //dz0holder=32.5;
    }
    else if(iPanel==3){//RR outer panel
        focal=1710.;// focal lenght
        xsiA1=1928.22;  //attaching points
        xsiA2=2652.77;
        etaA1=519.60;
        etaA2=1002.40;
        Dx=sqrt(pow(xsiA2-xsiA1,2.)+pow(etaA2-etaA1,2.))/2.;
        Dy=996./2.;
        //dz0holder=32.5;
    }
    else{//flat mirror or water surface
        zMirror=ui->dSB_zMirror->value();//z value (height) of the mirror/water surface
        Dx=300;
        Dy=200;
    }

    //j,i index of S cells containing the attaching points
    jiAP[0][0][0]=NINT(-Dy/DS)+jScenter;//index (j,i) of the cell at the attaching points
    jiAP[0][0][1]=NINT(-Dx/DS)+iScenter;//
    jiAP[0][0][2]=1;                       //               ^ X          ^ i
    jiAP[1][0][0]=NINT(Dy/DS)+jScenter;    //     [1][1][k] | [0][1][k]  |
    jiAP[1][0][1]=NINT(-Dx/DS)+iScenter;   //               |            |
    jiAP[1][0][2]=1;                       //       Y <-----0            |
    jiAP[0][1][0]=NINT(-Dy/DS)+jScenter;   //      (jScenter,iScenter)   |
    jiAP[0][1][1]=NINT(Dx/DS)+iScenter;    //                            |
    jiAP[0][1][2]=1;                       //     [1][0][k]   [0][0][k]  |
    jiAP[1][1][0]=NINT(Dy/DS)+jScenter;//   j <-----------------------
    jiAP[1][1][1]=NINT(Dx/DS)+iScenter;//        k=0->j     k=1->i
    jiAP[1][1][2]=1;                       //  i  CAM #1             CAM #2
    alpha=atan((etaA2-etaA1)/(xsiA2-xsiA1));//rotation angle to overlap xsi,eta -> x,z
    ui->dSB_focal->setValue(focal);
    //focal=ui->dSB_focal->value();
    //ui->doubleSpinBox_dz0hold->setValue(dz0holder);
    dz0holder=ui->doubleSpinBox_dz0hold->value();//NB per RondaReflex 95.20
    //double dzMirHold=ui->doubleSpinBox_dMirHold->value();
    //double zDeep=dz0holder-dzMirHold;
    //printf("Dx=%f zDeep=%f\n",Dx,zDeep);
    double xsi0,eta0;
    if(iPanel<2){
        xsi0=(xsiA1+xsiA2)/2.;
        eta0=(etaA1+etaA2)/2.;
    }
    else if(iPanel<4){
        double originLabRF=460.;//mm distance of LabRF origin form the attaching points closer to the parabola vertex
        xsi0=xsiA1+originLabRF*cos(alpha);//centro scacchiera collocato a distanza fissa da x attacchi lato vertice
        eta0=etaA1+originLabRF*sin(alpha);
        jiAP[0][0][1]=NINT(-originLabRF/DS)+iScenter;
        jiAP[1][0][1]=NINT(-originLabRF/DS)+iScenter;
        jiAP[0][1][1]=NINT((-originLabRF+979.2)/DS)+iScenter;
        jiAP[1][1][1]=NINT((-originLabRF+979.2)/DS)+iScenter;
    }
    if(iPanel<4){
        printf("xsi0=%f eta0=%f\n",xsi0,eta0);
        double xsiL=xsi0-dz0holder*sin(alpha);
        double etaL=eta0+dz0holder*cos(alpha);
        double xLab,xsi1,eta1;
        printf("xsiL=%f etaL=%f alpha=%f\n\tfocal=%f\n",xsiL,etaL,alpha/acos(-1)*180.,focal);
        if(iWr>0)
            printf("i xLab xsi eta xP zP tanXZ\n");
        double A,B,C,xsi,eta,tanXZ,xP,zP;
        for(int i=0; i<Ni;i++){
            xLab=double(i-iScenter)*DS;//in x,z LAB
            xsi1=xLab*cos(alpha)+xsiL;
            eta1=xLab*sin(alpha)+etaL;
            A=0.25/focal;
            B=1./tan(alpha);
            C=-eta1-xsi1/tan(alpha);
            xsi=(-B+sqrt(B*B-4.*A*C))/(2.*A);
            eta=0.25/focal*xsi*xsi;
            xP=(xsi-xsiL)*cos(alpha)+(eta-etaL)*sin(alpha);
            zP=-(xsi-xsiL)*sin(alpha)+(eta-etaL)*cos(alpha);
            tanXZ=tan(atan(0.5*xsi/focal)-alpha);
            if(iWr>0)
                printf("%d %f %f %f %f %f %f\n",i,xLab,xsi,eta,xP,zP,tanXZ);
            for(int j=0; j<Nj; j++){
                S[j][i][3]=zP;//ideal Z
                S[j][i][7]=zP;//experimental Z
                S[j][i][5]=tanXZ;//slopeX exp
                S[j][i][8]=tanXZ;//slopeX
                S[j][i][6]=0.;  //slopeY exp
                S[j][i][9]=0.;  //slopeY
            }
        }
        //focal position and unit solar versor in Lab Frame
        xsi=0.;
        eta=focal;
        pF[0]=(xsi-xsiL)*cos(alpha)+(eta-etaL)*sin(alpha);//focal line position
        pF[1]=0.;
        pF[2]=-(xsi-xsiL)*sin(alpha)+(eta-etaL)*cos(alpha);
        printf("focal position in Lab: pF=(%f ,%f , %f)\n",pF[0],pF[1],pF[2]);
        printf("idealHight completed!\n");
        fflush(stdout);
    }
    else{
        jiAP[0][0][0]=jScenter;
        jiAP[0][0][1]=iScenter;
        jiAP[1][0][0]=0;
        jiAP[1][0][1]=0;
        jiAP[1][0][2]=0;
        jiAP[0][1][0]=0;
        jiAP[0][1][1]=0;
        jiAP[0][1][2]=0;
        jiAP[1][1][0]=0;
        jiAP[1][1][1]=0;
        jiAP[1][1][2]=0;
        for(int i=0; i<Ni;i++){
            for(int j=0; j<Nj; j++){
                S[j][i][3]=zMirror;//ideal Z
                S[j][i][7]=zMirror;//experimental Z
                S[j][i][5]=0.;//slopeX exp
                S[j][i][8]=0.;//slopeX
                S[j][i][6]=0.;//slopeY exp
                S[j][i][9]=0.;//slopeY
            }
        }
    }
}


void kernelpro(int icam){
    printf("kernelpro icam=%d\n",icam);
    fflush(stdout);
    char filef[200];
    double DlabCam;
    Mat img,imgc,im_with_keypoints;
    if(icam==1)
        DlabCam=Pc_1[0]-x4set;
    else
        DlabCam=Pc_2[0]-x4set;
    int imgExist=1;
    int iFrame=-1;
    while(imgExist==1){//for(int iFrame=0; iFrame<nFrame; iFrame++){
        iFrame++;
        if(iWr==0){
            printf("\r\tiFrame=%d",iFrame);
            fflush(stdout);
        }
        else
        printf("Cam%d iFrame=%d x=%f xLabCam=%f\n",
               icam,iFrame,Xstart+Xstep*iFrame,Xstart+Xstep*iFrame+DlabCam);

        //            file=dirMeas+"/imgCam"+QString::number(icam)+"_"+QString::number(iFrame+1)+".bmp";
        //            img=imread(file.toStdString());

        //            file=dirMeas+"/"+QString::number(icam-1)+"/imgCam"+QString::number(icam)+"_"+QString::number(iFrame)+".bmp";
        //            img=imread(file.toStdString());
        QString file=dirMeas+"/"+QString::number(icam-1)+"/live-";
        sprintf(filef,"%s%010d.bmp",(file.toStdString()).c_str(),iFrame);
        QFileInfo check_file(filef);
        // check if file exists and if yes: Is it really a file and no directory?
        if (check_file.exists() && check_file.isFile()) {
            img=imread(filef);
        } else {
            imgExist=0;
            continue;
        }

        if(icam==1){
            if(undst1==1)
                undistort(img,imgc,cameraMatrix_1,distCoeffs_1);
            else
                imgc=img;
        }
        else{
            if(undst2==1)
                undistort(img,imgc,cameraMatrix_2,distCoeffs_2);
            else
                imgc=img;
        }

        // set up and create the detector using the parameters
        cv::Ptr<cv::SimpleBlobDetector> blob_detector = cv::SimpleBlobDetector::create(params);

        // detect!
        vector<cv::KeyPoint> keypoints;
        blob_detector -> detect(imgc, keypoints);

        if(iWr==1 || iWr==2){
            // Draw detected blobs as red circles.
            // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures the size of the circle corresponds to the size of blob
            drawKeypoints( imgc, keypoints, im_with_keypoints, Scalar(0,0,255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );
            if(icam==1)
                imshow("Camera1",im_with_keypoints);
            else
                imshow("Camera2",im_with_keypoints);
            //imwrite(dir.toStdString()+"/imageWithBlob.JPEG",im_with_keypoints);
        }

        // extract the x y coordinates of the keypoints:
        int nBlobs = keypoints.size();
        //printf("Found  nBlobs = %d\n",nBlobs);
        if(nBlobs>0){
            int iw=0;
            double blobRes[nBlobs][3];
            for (int i=0; i<nBlobs; i++){
                float J = keypoints[i].pt.x;
                float I = keypoints[i].pt.y;
                float S = keypoints[i].size;
                iw=i;
                if(i>0){
                    while(J < blobRes[iw-1][0]){
                        iw--;
                    }
                    for(int ii=i-1;ii>=iw;ii--){
                        blobRes[ii+1][0]=blobRes[ii][0];
                        blobRes[ii+1][1]=blobRes[ii][1];
                        blobRes[ii+1][2]=blobRes[ii][2];
                    }
                }
                blobRes[iw][0]=J;
                blobRes[iw][1]=I;
                blobRes[iw][2]=S;
            }
            int ibC100=-1,iB=0;
            double Dj,Di,Djsum=0.,DjMax=-1000.,DjMin=1000.;
            if(iWr==1 || iWr==2)
            printf("Found nBlobs=%d ",nBlobs);
            for (int ib=0; ib<nBlobs-1; ib++){
                Dj=blobRes[iB+1][0]-blobRes[iB][0];
                Di=blobRes[iB+1][1]-blobRes[iB][1];
                if(sqrt(Dj*Dj+Di*Di)<minDis){//merge the two blobs
                    blobRes[iB][0]=0.5*(blobRes[iB][0]+blobRes[iB+1][0]);
                    blobRes[iB][1]=0.5*(blobRes[iB][1]+blobRes[iB+1][1]);
                    for(int k=iB+1;k<nBlobs-1;k++){//move 1 step back
                        for(int kk=0;kk<3;kk++)
                            blobRes[k][kk]=blobRes[k+1][kk];
                    }
                    nBlobs--;
                    if(iWr==1 || iWr==2)
                    printf("two blobs were merged: %d with %d\nNew nBlobs=%d\n",iB,iB+1,nBlobs);
                    continue;
                }
                //printf("iBlob = %d j=%f i=%f S=%f Di= %f Dj= %f\n",
                //       iB+1,blobRes[iB+1][0],blobRes[iB+1][1],blobRes[iB+1][2],Di,Dj);
                if(ibC100==-1 && Di<-5. && (blobRes[iB+2][1]-blobRes[iB+1][1])>10){
                    ibC100=iB+1;//iCenter== point source N. 100
                    if(iWr==1 || iWr==2)
                    printf("found ibC100=%d ",ibC100);
                }
                else if(ibC100==-1 && Di>5. && (blobRes[iB+2][1]-blobRes[iB+1][1])<-10){
                    double c100[3];
                    for(int k=0;k<3;k++)
                        c100[k]=blobRes[iB+2][k];
                    for(int k=0;k<3;k++)
                        blobRes[iB+2][k]=blobRes[iB+1][k];
                    for(int k=0;k<3;k++)
                        blobRes[iB+1][k]=c100[k];
                    ibC100=iB+1;//iCenter== point source N. 100
                    if(iWr==1 || iWr==2){
                        printf("found ibC100=%d ",ibC100);
                        printf("\tiBlob = %d j=%f i=%f\n",
                               iB+1,blobRes[iB+1][0],blobRes[iB+1][1]);
                        printf("\tiBlob = %d j=%f i=%f\n",
                               iB+2,blobRes[iB+2][0],blobRes[iB+2][1]);
                    }
                }
                else if(abs(Di)<5.){
                    DjMin=min(DjMin,Dj);
                    DjMax=max(DjMax,Dj);
                    Djsum=Djsum+Dj;
                }
                iB++;
            }
            Dj=Djsum/double(nBlobs-2);
            //                    if(DjMax>1.75*Dj || DjMin<0.4*Dj)
            //                        iW=1;
            //                    else
            //                        iW=0;
            if(iWr==1 || iWr==2)
            printf("<Dj>=%f DjMin=%f DjMax=%f\n",
                   Dj,DjMin,DjMax);
            if(ibC100!=-1){
                int nSource;
                int ib;
                int iStep=-1;
                for (int i=0; i<nBlobs; i++){
                    ib=ibC100-i;
                    if(ib<0){
                        ib=ibC100+(i-ibC100);
                        iStep=1;
                    }
                    if(ib == ibC100 || ib == ibC100+1){
                        if(blobRes[ib][1] > blobRes[ibC100][1])
                            nSource=100;
                        else
                            nSource=99;
                    }
                    //    printf("i=%d ib=%d ",i,ib);
                    iTot++;
                    if(ib-iStep>=0 && ib-iStep<nBlobs && abs(blobRes[ib][0]-blobRes[ib-iStep][0]) > 1.75*abs(Dj)){
                        int DnS=NINT(abs((blobRes[ib][0]-blobRes[ib-iStep][0])/Dj));
                        nSource=nSource+DnS*double(-iStep);
                    }
                    else
                        nSource=nSource-iStep;

                    results[iTot*6]   = iFrame; //N frame
                    results[iTot*6+1] = icam;   //N camera
                    results[iTot*6+2] = nSource;//N source
                    results[iTot*6+4] = blobRes[ib][0];//j bubble-coordinate (px)
                    results[iTot*6+5] = blobRes[ib][1];//i bubble-coordinate (px)
                    if(iWr==2){
                        double ddjj;
                        if(ib-iStep>=0 && ib-iStep<nBlobs)
                            ddjj=blobRes[ib][0]-blobRes[ib-iStep][0];
                        else
                            ddjj=69;
                        printf(" nS=%d j=%f i=%f S=%f Dj=%f\n",
                               nSource,blobRes[ib][0],blobRes[ib][1],
                                blobRes[ib][2],ddjj);
                        if(abs(abs(ddjj)-abs(DjMin))<0.001)
                            printf("^^^^^^^^^^^^ MIN  !!!!!!!!!\n");
                        if(abs(abs(ddjj)-abs(DjMax))<0.001)
                            printf("^^^^^^^^^^^^ MAX  !!!!!!!!!\n");
                        fflush(stdout);
                    }
                }
                if(iWr==1 || iWr==2)
                    waitKey(0);
            }
            else{
                printf("ATTENTION: invalid iCenter!\nprocess aborted!\n"
                       "Press return to continue\n");
                waitKey(0);
                return;
            }
        }
    }
    printf("\tkernelpro completed with iTot=%d\n",iTot);
    fflush(stdout);
}

void slopeComputing(){
    printf("Slope computing on %d data .....",iTot);
    fflush(stdout);
    int jr,ir,j,i,nSource,iCam,iFrame;
    double t,Dt;
    double x0,y0,z0,vr[3],x,y,z,xS,yS,zS;
    for(j=0;j<Nj;j++)
        for(i=0;i<Ni;i++)
            for(int k=0;k<5;k++){
                if(k==3) k++;
                S[j][i][k]=0.;
            }
    for(int idato=0;idato<=iTot;idato++){
        //printf("\ridato=%d",idato);
        //fflush(stdout);
        iFrame=NINT(results[idato*6]);
        iCam=NINT(results[idato*6+1]);
        nSource=NINT(results[idato*6+2]);
        double dZcam=0.,dPitch=0.,dRoll=0.;
        if(iWr>0)
            printf("idato=%d -> iFrame=%d iCam=%d nSource=%d\n",idato,iFrame,iCam,nSource);
        if(iArtCor==1){
            int iLastArtCor=215;
            double xMotor=Xstart+Xstep*iFrame;
            while(corArte[iLastArtCor][0]>xMotor && iLastArtCor>0){
                iLastArtCor--;
            }
            dZcam=corArte[iLastArtCor][iCam];
            dPitch=corArte[iLastArtCor][3];
            dRoll=corArte[iLastArtCor][4];
            if(iWr==1 || iWr==2){
                printf("xMotor=%f dZcam=%f dPitch=%f dRoll=%f\n",
                       corArte[iLastArtCor][0],dZcam,dPitch,dRoll);
            }
        }
        double fx,fy,cx,cy,PA[6],DlabCam;
        if(iCam==1){
            cx=cx_1;
            cy=cy_1;
            fx=fx_1;
            fy=fy_1;
            for(int i=0;i<3;i++){
                PA[i]=Pc_1[i];
                PA[3+i]=EA_1[i];
            }
            DlabCam=Pc_1[0]-x4set;
        }
        else{
            cx=cx_2;
            cy=cy_2;
            fx=fx_2;
            fy=fy_2;
            for(int i=0;i<3;i++){
                PA[i]=Pc_2[i];
                PA[3+i]=EA_2[i];
            }
            DlabCam=Pc_2[0]-x4set;
        }
        pxX=results[idato*6+4]-cx; //j bubble-coordinate (px)
        pxY=results[idato*6+5]-cy; //i bubble-coordinate (px)
        unitV12(0.0,0.0,0.0,pxX/fx,pxY/fy,1.);
        if(iWr==2){
            printf("PA: %f %f %f %f %f %f\n",PA[0],PA[1],PA[2],PA[3],PA[4],PA[5]);
            printf("Vr RifCam:\t%f\t%f\t%f\n",Vservice[0],Vservice[1],Vservice[2]);
        }
        Plane2Earth(PA[3],PA[4],PA[5],Vservice[0],Vservice[1],Vservice[2]);//trasformazione di uvr da rifDrone a rifLab
        if(iWr==2)
            printf("Vr RifLab:\t%f\t%f\t%f\n",Vservice[0],Vservice[1],Vservice[2]);
        if(iArtCor==1)
            Plane2Earth(0.,-dPitch,-dRoll,Vservice[0],Vservice[1],Vservice[2]);
        for(int k=0;k<3;k++)
            vr[k]=Vservice[k];
        x0=Xstart+Xstep*iFrame+DlabCam;//camera x coordinate
        y0=PA[1];
        z0=PA[2]+dZcam;
        if(iWr>0)
            printf("idato=%d nSource=%d Pcam(%f, %f, %f) Vr(%f, %f, %f)\n",
                idato,nSource,x0,y0,z0,vr[0],vr[1],vr[2]);
        if(vr[0]!=0. || vr[1]!=0.){
            if(abs(vr[0])>abs(vr[1]))
                Dt=DS/abs(vr[0]);//step needed to change cell
            else
                Dt=DS/abs(vr[1]);//step needed to change cell
            t=abs(z0/vr[2])*0.8;
            x=vr[0]*t+x0;//point P on the reflective surface
            y=vr[1]*t+y0;
            z=vr[2]*t+z0;
            i=NINT(x/DS)+iScenter;
            j=NINT(y/DS)+jScenter;
            //printf("t=%f Dt=%f\n",t,Dt);
            double DelZold,DelZ;
            t=t-Dt;
            do{
                t=t+Dt;
                DelZold=abs(z-S[j][i][7]);
                x=vr[0]*(t+Dt)+x0;//point P on the reflective surface
                y=vr[1]*(t+Dt)+y0;
                z=vr[2]*(t+Dt)+z0;
                i=NINT(x/DS)+iScenter;
                j=NINT(y/DS)+jScenter;
                DelZ=z-S[j][i][7];
                //printf("\tt=%f Dt=%f iS=%d jS=%d DelZ=%f DelZold=%f\n",
                //       t+Dt,Dt,i,j,DelZ,DelZold);
            }while(DelZ>0.);//while(j>=0 && j <Nj && i>=0 && i<Ni && z>S[j][i][7]);
            t=t+Dt*DelZold/(DelZold-DelZ);
        }
        else{
            i=NINT(x0/DS)+iScenter;
            j=NINT(y0/DS)+jScenter;
            t=z0-S[j][i][7];
        }
        if(j>=0 && j <Nj && i>=0 && i<Ni ){
            x=vr[0]*t+x0;//point P on the reflective surface
            y=vr[1]*t+y0;
            z=vr[2]*t+z0;
            ir=NINT(x/DS)+iScenter;
            jr=NINT(y/DS)+jScenter;
            //printf("finale: t=%f iS=%d jS=%d MF=%f\n",t,i,j,z-S[jr][ir][7]);
            Plane2Earth(Es[0],Es[1],Es[2],source[nSource][0],source[nSource][1],source[nSource][2]);
            xS=Vservice[0]+Ps[0];
            yS=Vservice[1]+Ps[1];
            zS=Vservice[2]+Ps[2];
            if(jr==30 && iWr==3){
                printf("\nidato= %d obtained by Cam #%d at frame #%d\n",idato,NINT(results[idato*6+1]),NINT(results[idato*6]));
                printf("Cam: %f %f %f\n",x0,y0,z0);
                printf("  P: %f %f %f\n",x,y,z);
                printf("  S: %f %f %f nSource= %d\n",xS,yS,zS,nSource);
                printf("Vr: %f %f %f\n",-vr[0],-vr[1],-vr[2]);
            }

            unitV12(xS,yS,zS,x,y,z);//unit vector source->its image
            if(jr==30 && iWr==3)
                printf("Vi: %f %f %f\n",Vservice[0],Vservice[1],Vservice[2]);
            sumV12(-Vservice[0],-Vservice[1],-Vservice[2],-vr[0],-vr[1],-vr[2]);//normal unit vector is returned in Vservice
            if(iWr==2){
                printf("Vn: %f %f %f in [%d][%d]\n",Vservice[0],Vservice[1],Vservice[2],jr,ir);
                printf("dz/dx= %f dz/dy=%f\n",-Vservice[0]/Vservice[2],-Vservice[1]/Vservice[2]);
                fflush(stdout);
            }
            if(jr>=0 && jr <Nj && ir>=0 && ir<Ni){
                if(int(S[jr][ir][4])==0){
                    S[jr][ir][0]=Vservice[0];
                    S[jr][ir][1]=Vservice[1];
                    S[jr][ir][2]=Vservice[2];
                }else{
                    S[jr][ir][0]=S[jr][ir][0]+Vservice[0];//the average value is computed in void inpainting()
                    S[jr][ir][1]=S[jr][ir][1]+Vservice[1];
                    S[jr][ir][2]=S[jr][ir][2]+Vservice[2];
                }
                S[jr][ir][4]++;
            }
        }
//        waitKey(0);
    }
    printf("\tslopeComputing done!\n");
    fflush(stdout);
}

void VISproPT::intFat(){
    printf("computing intFat...");
    double phiS=0.0047;//(rad) solar divergence half-angle
    double Rrec=35.;//(mm) receiver radius
    double thetaL=ui->dSB_thetaL->value();
    thetaL=thetaL/180.*pig;
    double uvN[3];//normal unit vector
    double uvR[3];//reflected unit vector
    double pP[3];//point of reflection
    //uvS in parabola natural frame
    uvS[0]=0.;
    uvS[1]=-sin(thetaL);
    uvS[2]=-cos(thetaL);
    Plane2Earth(0.,-alpha,0.,uvS[0],uvS[1],uvS[2]);
    for(int k=0;k<3;k++)
        uvS[k]=Vservice[k];//uvS in Lab frame
    for(int i=isMIN;i<=isMAX;i++){
        for(int j=jsMIN;j<=jsMAX;j++){
            S[j][i][10]=0.;
            if(NINT(S[j][i][4])>0){
                double psca=0.;
                double mod=0.;
                for(int k=0; k<3; k++){
                    uvN[k]=S[j][i][k];
                    psca=psca+uvN[k]*uvS[k];//scalar multiplication
                    mod=mod+uvN[k]*uvN[k];
                }
                pP[0]=(i-iScenter)*DS;
                pP[1]=0;//(j-jScenter)*DS; y value is not relevant
                pP[2]=S[j][i][7];
                sumV12(uvS[0],uvS[1],uvS[2],-2.*psca*uvN[0],-2.*psca*uvN[1],-2.*psca*uvN[2]);//reflected unit vector
                for(int k=0; k<3; k++)
                    uvR[k]=Vservice[k];
//                if(j==75){
//                    printf("j=%d\ti=%d NdatCell=%d\n",j,i,NINT(S[j][i][4]));
//                    printf("uvS=(%f , %f , %f)\n",uvS[0],uvS[1],uvS[2]);
//                    printf("uvN=(%f , %f , %f)\n",uvN[0],uvN[1],uvN[2]);
//                    printf("uvR=(%f , %f , %f)\n",uvR[0],uvR[1],uvR[2]);
//                }
                //distance reflected ray from focus line
                double a=uvR[0]*uvR[0]+uvR[2]*uvR[2];
                double b=2.*(uvR[0]*(pP[0]-pF[0])+uvR[2]*(pP[2]-pF[2]));
                double c=pow((pP[0]-pF[0]),2.)+pow((pP[2]-pF[2]),2.);
                //double delta=b*b-4.*a*c;
                double d=sqrt(-b*b/4./a+c);// distance by setting delta=0
                double t=-b/(2.*a); //path length

                double Rbeam=t*phiS;//radius of reflected solar beam
                double d1=d-Rbeam;
                double d2=d+Rbeam;
                double theta,IntFactor,lost;
//                if(j==75)
//                    printf("pP=(%f , %f ,%f) t=%f d=%f d1=%f d2=%f\n",
//                       pP[0],pP[1],pP[2],t,d,d1,d2);
                if(d1>=Rrec)
                    IntFactor=0.;
                else{
                    lost=0.;
                    if(d1<-Rrec){
                        theta=acos((d+Rrec)/Rbeam);
                        lost=(theta-sin(theta)*cos(theta));
                    }
                    if(d2>Rrec){
                        theta=acos((Rrec-d)/Rbeam);
                        lost=lost+(theta-sin(theta)*cos(theta));
                    }
                    IntFactor=1.-lost/pig;
                }

                S[j][i][10]=IntFactor;
//                if(j==75){
//                    printf("IntFact=%f \n",IntFactor);
//                    //waitKey(0);
//                }
            }
        }
    }
}


int kernelSetCam(int iSetCamMeth){
    printf("kernelSetCam with Nt=%d NbCam1=%d NbCam2=%d zMirror=%f...\n",
           Ncorn,NbCam1,NbCam2,zMirror);
    fflush(stdout);
    //specular values
    Ps[2]=-Ps[2];//only z
    for(int i=1;i<3;i++)//except yaw
        Es[i]= -Es[i];
    int info;
    int m=2*(Ncorn+Ncorn)+2*(NbCam1+NbCam2);
    int n;
    if(iSetCamMeth==0)
        n=16;
    else
        n=10;
    int lwa=m*n+5*n+m;
    //    int iwa[n];
    //    double xx[n],fvec[m],wa[lwa];
    int* iwa=nullptr;
    iwa = new int[n];
    double* xx=nullptr;
    xx = new double[n];
    double* fvec=nullptr;
    fvec = new double[m];
    double* wa=nullptr;
    wa = new double[lwa];
    double tol=sqrt(dpmpar(1));
    //pTF2[0].Nt=Nt;
    xx[0]=Pc_1[0];
    xx[1]=Pc_1[1];
    xx[2]=Pc_1[2];
    xx[3]=EA_1[0];
    xx[4]=EA_1[1];
    xx[5]=EA_1[2];
    //xx[6]=Pc_2[0];
    xx[6]=Pc_2[1];
    //xx[8]=Pc_2[2];
    xx[7]=EA_2[0];
    xx[8]=EA_2[1];
    xx[9]=EA_2[2];
    if(n==16){
        xx[10]=fx_1;
        xx[11]=fx_2;
        xx[12]=cx_1;
        xx[13]=cy_1;
        xx[14]=cx_2;
        xx[15]=cy_2;
    }
//    xx[18]=Es[0];
//    xx[ ]=Es[1];
//    xx[19]=Es[2];
//    xx[20]=distCoeffs_1.at<double>(0);
//    xx[21]=distCoeffs_1.at<double>(1);
//    xx[22]=distCoeffs_1.at<double>(2);
//    xx[23]=distCoeffs_1.at<double>(3);
//    xx[24]=distCoeffs_1.at<double>(4);
//    xx[25]=distCoeffs_2.at<double>(0);
//    xx[26]=distCoeffs_2.at<double>(1);
//    xx[27]=distCoeffs_2.at<double>(2);
//    xx[28]=distCoeffs_2.at<double>(3);
//    xx[29]=distCoeffs_2.at<double>(4);
    printf("\tbest-fit: n=%d parameters and m=%d data Ncorn=%d\n",n,m,Ncorn);

    info=lmdif1(fcn, &pTF2, m, n, xx,fvec, tol, iwa, wa, lwa);

    Pc_1[0]=xx[0];
    Pc_1[1]=xx[1];
    Pc_1[2]=xx[2];
    EA_1[0]=xx[3];
    EA_1[1]=xx[4];
    EA_1[2]=xx[5];
    //Pc_2[0]=xx[6];
    Pc_2[0]=Pc_1[0]+Dx21;
    Pc_2[1]=xx[6];
    //Pc_2[2]=xx[8];
    Pc_2[2]=Pc_1[2]+Dz21;
    EA_2[0]=xx[7];
    EA_2[1]=xx[8];
    EA_2[2]=xx[9];
    if(n==16){
        fx_1=xx[10];
        fx_2=xx[11];
        cx_1=xx[12];
        cy_1=xx[13];
        cx_2=xx[14];
        cy_2=xx[15];
    }
//    Es[0]=xx[18];
//    Es[1]=xx[4];
//    Es[2]=xx[19];
//    distCoeffs_1.at<double>(0)=xx[20];
//    distCoeffs_1.at<double>(1)=xx[21];
//    distCoeffs_1.at<double>(2)=xx[22];
//    distCoeffs_1.at<double>(3)=xx[23];
//    distCoeffs_1.at<double>(4)=xx[24];
//    distCoeffs_2.at<double>(0)=xx[25];
//    distCoeffs_2.at<double>(1)=xx[26];
//    distCoeffs_2.at<double>(2)=xx[27];
//    distCoeffs_2.at<double>(3)=xx[28];
//    distCoeffs_2.at<double>(4)=xx[29];
    Ps[2]=-Ps[2];//restore original value (only z)
    for(int i=1;i<3;i++)//except yaw
        Es[i]= -Es[i];
    if(wa){
        delete[] wa;
        wa = nullptr;
    }
    if(fvec){
        delete [] fvec;
        fvec=nullptr;
    }
    if(xx){
        delete [] xx;
        xx=nullptr;
    }
    if(iwa){
        delete []  iwa;
        iwa=nullptr;
    }
    printf("\tPcam1= %f\t%f\t%f\tEcam1= %f\t%f\t%f\n",
           Pc_1[0],Pc_1[1],Pc_1[2],EA_1[0]*180./pig,EA_1[1]*180./pig,EA_1[2]*180./pig);
    printf("\tPcam2= %f\t%f\t%f\tEcam1= %f\t%f\t%f\n",
           Pc_2[0],Pc_2[1],Pc_2[2],EA_2[0]*180./pig,EA_2[1]*180./pig,EA_2[2]*180./pig);
    return(info);
}


int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag){
    /* calculate the functions at x and return the values in fvec[0] through fvec[m-1] */
    struct pointToFit2 *pTF2 = (struct pointToFit2 *)p;
    int fredeg=m-n;
    int Nt=Ncorn;
    int k=0;
    double chi2=0.;
    Pc_1[0]=x[0];
    Pc_1[1]=x[1];
    Pc_1[2]=x[2];
    EA_1[0]=x[3];
    EA_1[1]=x[4];
    EA_1[2]=x[5];
    //Pc_2[0]=x[6];
    Pc_2[0]=Pc_1[0]+Dx21;
    Pc_2[1]=x[6];
    //Pc_2[2]=x[8];
    Pc_2[2]=Pc_1[2]+Dz21;
    EA_2[0]=x[7];
    EA_2[1]=x[8];
    EA_2[2]=x[9];
    if(n==16){
        fx_1=x[10];
        fx_2=x[11];
        cx_1=x[12];
        cy_1=x[13];
        cx_2=x[14];
        cy_2=x[15];
    }
//    Es[0]=x[18];
//    Es[1]=xx[4];
//    Es[2]=x[19];
//    distCoeffs_1.at<double>(0)=x[20];
//    distCoeffs_1.at<double>(1)=x[21];
//    distCoeffs_1.at<double>(2)=x[22];
//    distCoeffs_1.at<double>(3)=x[23];
//    distCoeffs_1.at<double>(4)=x[24];
//    distCoeffs_2.at<double>(0)=x[25];
//    distCoeffs_2.at<double>(1)=x[26];
//    distCoeffs_2.at<double>(2)=x[27];
//    distCoeffs_2.at<double>(3)=x[28];
//    distCoeffs_2.at<double>(4)=x[29];

    //chessboard data
    for(int i=0;i<Nt;i++){
        if(px[i][0] < 0.9e+9){
            xyz2pxpy(1,P[i][0],P[i][1],P[i][2]);//compute pxX pxY given X,Y,Z,Pc[3],EA[3]
            //PxPyDistorted(1,pxX,pxY);
            px[i][2]=pxX+cx_1;
            px[i][3]=pxY+cy_1;
            fvec[k]=cx_1+pxX-px[i][0];
            chi2=chi2+fvec[k]*fvec[k];
            k++;
            fvec[k]=cy_1+pxY-px[i][1];
            chi2=chi2+fvec[k]*fvec[k];
            k++;
        }
    }
    for(int i=Nt;i<2*Nt;i++){
        if(px[i][0] < 0.9e+9){
            xyz2pxpy(2,P[i][0],P[i][1],P[i][2]);//compute pxX pxY given X,Y,Z,Pc[3],EA[3]
            //PxPyDistorted(2,pxX,pxY);
            px[i][2]=pxX+cx_2;
            px[i][3]=pxY+cy_2;
            fvec[k]=cx_2+pxX-px[i][0];
            chi2=chi2+fvec[k]*fvec[k];
            k++;
            fvec[k]=cy_2+pxY-px[i][1];
            chi2=chi2+fvec[k]*fvec[k];
            k++;
        }
    }

    //source reflected on flat mirror thick=zMirror
    Ps[2]=Ps[2]+zMirror;//correction for the Mirror thickness
    int nSource=0;
    for(int i=0;i<NbCam1;i++){
        nSource=NINT(Bcam1[i][2]);
        //trasformation point source rifSource -> rifLab
        Plane2Earth(Es[0],Es[1],Es[2],
                source[nSource][0],source[nSource][1],source[nSource][2]);
        for(int ii=0;ii<3;ii++)
            Vservice[ii]=Vservice[ii]+Ps[ii];
        //printf("iCam1=%d x,y,z= %f\t%f\t%f\n",i,Vservice[0],Vservice[1],Vservice[2]);
        xyz2pxpy(1,Vservice[0],Vservice[1],Vservice[2]);//compute pxX pxY given X,Y,Z,Pc[3],EA[3]
        //PxPyDistorted(1,pxX,pxY);
        //px[2*Nt+i][0]=Bcam1[i][0];
        //px[2*Nt+i][1]=Bcam1[i][1];
        px[2*Nt+i][2]=pxX+cx_1;
        px[2*Nt+i][3]=pxY+cy_1;
        //printf("pxX=%f <-> %f\tpxY=%f <->%f\n",pxX,Bcam1[i][0],pxY,Bcam1[i][1]);
        fvec[k]=cx_1+pxX-Bcam1[i][0];
        chi2=chi2+fvec[k]*fvec[k];
        k++;
        fvec[k]=cy_1+pxY-Bcam1[i][1];
        chi2=chi2+fvec[k]*fvec[k];
        k++;
    }
    for(int i=0;i<NbCam2;i++){
        nSource=NINT(Bcam2[i][2]);
        //trasformation point source rifSource -> rifLab
        Plane2Earth(Es[0],Es[1],Es[2],
                source[nSource][0],source[nSource][1],source[nSource][2]);
        for(int ii=0;ii<3;ii++)
            Vservice[ii]=Vservice[ii]+Ps[ii];
        //printf("iCam2=%d x,y,z= %f\t%f\t%f\n",i,Vservice[0],Vservice[1],Vservice[2]);
        xyz2pxpy(2,Vservice[0],Vservice[1],Vservice[2]);//compute pxX pxY given X,Y,Z,Pc[3],EA[3]
        //PxPyDistorted(2,pxX,pxY);
        //px[2*Nt+NbCam1+i][0]=Bcam2[i][0];
        //px[2*Nt+NbCam1+i][1]=Bcam2[i][1];
        px[2*Nt+NbCam1+i][2]=pxX+cx_2;
        px[2*Nt+NbCam1+i][3]=pxY+cy_2;
        //printf("pxX=%f <-> %f\tpxY=%f <->%f\n",pxX,Bcam2[i][0],pxY,Bcam2[i][1]);
        fvec[k]=cx_2+pxX-Bcam2[i][0];
        chi2=chi2+fvec[k]*fvec[k];
        k++;
        fvec[k]=cy_2+pxY-Bcam2[i][1];
        chi2=chi2+fvec[k]*fvec[k];
        k++;
    }
    Ps[2]=Ps[2]-zMirror;//remove correction
    chi2fin=chi2/double(fredeg);
    return(0);
}


int fcnRefCam(void *p, int m, int n, const double *x, double *fvec, int iflag){
    /* calculate the functions at x and return the values in fvec[0] through fvec[m-1] */
    struct pointToFit2 *pTF2 = (struct pointToFit2 *)p;
    int iRefOpt=NINT(pTF2[0].Nt);
    if(iRefOpt==0){
        fx_1=x[0];
        fy_1=x[0];
        cx_1=x[1];
        cy_1=x[2];
        fx_2=x[3];
        fy_2=x[3];
        cx_2=x[4];
        cy_2=x[5];
    }
    else if(iRefOpt==1){
        Pc_1[0]=x[0];
        Pc_1[1]=x[1];
        Pc_1[2]=x[2];
        EA_1[0]=x[3];
        EA_1[1]=x[4];
        EA_1[2]=x[5];
        //Pc_2[0]=x[6];
        Pc_2[0]=Pc_1[0]+Dx21;
        Pc_2[1]=x[6];
        //Pc_2[2]=x[8];
        Pc_2[2]=Pc_1[2]+Dz21;
        EA_2[0]=x[7];
        EA_2[1]=x[8];
        EA_2[2]=x[9];
    }
    else{
        fx_1=x[0];
        fy_1=x[0];
        cx_1=x[1];
        cy_1=x[2];
        fx_2=x[3];
        fy_2=x[3];
        cx_2=x[4];
        cy_2=x[5];
        Pc_1[0]=x[6];
        Pc_1[1]=x[7];
        Pc_1[2]=x[8];
        EA_1[0]=x[9];
        EA_1[1]=x[10];
        EA_1[2]=x[11];
        //Pc_2[0]=x[12];
        Pc_2[0]=Pc_1[0]+Dx21;
        Pc_2[1]=x[12];
        //Pc_2[2]=x[14];
        Pc_2[2]=Pc_1[2]+Dz21;
        EA_2[0]=x[13];
        EA_2[1]=x[14];
        EA_2[2]=x[15];
    }
    slopeComputing();//new slope computing
    //averaging
    for(int i=isMIN; i<=isMAX; i++){
        for(int j=jsMIN; j<=jsMAX; j++){
            if(NINT(S[j][i][4])>0){
                double mod=0.;
                for(int k=0;k<3;k++)//averaging
                    mod=mod+S[j][i][k]*S[j][i][k];
                for(int k=0;k<3;k++)
                    S[j][i][k]=S[j][i][k]/sqrt(mod);
                S[j][i][5]=-S[j][i][0]/S[j][i][2];//dz/dx (partial x-derivative)
                S[j][i][6]=-S[j][i][1]/S[j][i][2];//dz/dy (partial y-derivative)
            }
        }
    }
    //inpainting
    for(int i=isMIN; i<=isMAX; i++){
        for(int j=jsMIN; j<=jsMAX; j++){
            if(NINT(S[j][i][4])==0){
                if(j-1>=jsMIN && j+1<=jsMAX){
                    if(S[j-1][i][4]>0 && S[j+1][i][4]>0){
                        for(int k=0; k<3 ; k++)
                            S[j][i][k]=S[j][i][k]+(S[j-1][i][k]+S[j+1][i][k])/2.;
                        S[j][i][4]++;
                    }
                }
                if(j-1>=jsMIN && j+1<=jsMAX && i-1 >=isMIN && i+1<=isMAX){
                    if(S[j-1][i-1][4]>0 && S[j+1][i+1][4]>0){
                        for(int k=0; k<3 ; k++){
                            S[j][i][k]=S[j][i][k]+(S[j-1][i-1][k]+S[j+1][i+1][k])/2.;
                            S[j][i][k]=S[j][i][k]+(S[j-1][i+1][k]+S[j+1][i-1][k])/2.;
                        }
                        S[j][i][4]++;
                    }
                }
                if(i-1>=isMIN && i+1<=isMAX){
                    if(S[j][i-1][4]>0 && S[j][i+1][4]>0){
                        for(int k=0; k<3 ; k++)
                            S[j][i][k]=S[j][i][k]+(S[j][i-1][k]+S[j][i+1][k])/2.;
                        S[j][i][4]++;
                    }
                }
                if(NINT(S[j][i][4])>0){
                    double mod=0.;
                    for(int k=0;k<3;k++)//averaging
                        mod=mod+S[j][i][k]*S[j][i][k];
                    for(int k=0; k<3 ; k++)
                        S[j][i][k]=S[j][i][k]/sqrt(mod);
                    S[j][i][5]=-S[j][i][0]/S[j][i][2];//dz/dx (partial x-derivative)
                    S[j][i][6]=-S[j][i][1]/S[j][i][2];//dz/dy (partial y-derivative)
                }
            }
        }
    }
    //fvec
    int k=0;
    double chi2=0.;
    for(int i=isMIN+5; i<=isMAX-5; i++){
        for(int j=jsMIN+5; j<=jsMAX-5; j++){
            if(NINT(S[j][i][4])>0){
                fvec[k]=(S[j][i][5]-S[j][i][8])*1.E+03;
                chi2=chi2+fvec[k]*fvec[k];
                k++;
                fvec[k]=(S[j][i][6]-S[j][i][9])*1.E+03;
                chi2=chi2+fvec[k]*fvec[k];
                k++;
            }else{
                fvec[k]=0.;
                k++;
                fvec[k]=0.;
                k++;
            }
        }
    }
    printf(">>>>>>chi2=%f iflag=%d \n",chi2,iflag);
    if(iRefOpt==0 || iRefOpt==2)
        printf("\t\tf_1=%f cx_1=%f cy_1=%f\n\t\tf_2=%f cx_2=%f cy_2=%f\n",
           x[0],x[1],x[2],x[3],x[4],x[5]);
    if(iRefOpt>0){
        printf("\tPcam1= %f\t%f\t%f\tEcam1= %f\t%f\t%f\n",
               Pc_1[0],Pc_1[1],Pc_1[2],EA_1[0]*180./pig,EA_1[1]*180./pig,EA_1[2]*180./pig);
        printf("\tPcam2= %f\t%f\t%f\tEcam1= %f\t%f\t%f\n",
               Pc_2[0],Pc_2[1],Pc_2[2],EA_2[0]*180./pig,EA_2[1]*180./pig,EA_2[2]*180./pig);
    }
    chi2fin=chi2;
    QFile fLog(dir+"/fLog.txt");
    if (!fLog.open(QIODevice::WriteOnly | QIODevice::Append))
        return(1);
    QTextStream out(&fLog);
    out<<"iCam1&2"<<"\t"<<chi2;
    for(int i=0;i<n;i++)
        out<<"\t"<<x[i];
    out<<"\n";
    fLog.close();
    fflush(stdout);
    return(0);
}



void InfoFit(int info){
    if(info==0)
        printf("Fit: 0= improper input parameters\n");
    else if(info==1)
        printf("Fit: 1= relative error in the sum of squares is at most tol\n");
    else if(info==2)
        printf("Fit: 2= relative error between x and the solution is at most tol\n");
    else if(info==3)
        printf("Fit: 3= conditions for info = 1 and info = 2 both hold\n");
    else if(info==4)
        printf("Fit: 4= fvec is orthogonal to the columns of the jacobian to machine precision\n");
    else if(info==5)
        printf("Fit: 5= number of calls to fcn has reached or exceeded 200*(n+1)\n");
    else if(info==6)
        printf("Fit: 6= tol is too small. No further reduction in SUM_sqare is possible\n");
    else if(info==7)
        printf("Fit: 7= tol is too small. no further improvement in the approximate solution x is possible\n");
    fflush(stdout);
}


//void pxpy2xyz(double di,double dj){ //compute XYZ in Lab Frame for a given pixel (i,j) on the laser trace
//   double vd[3],t;
//   vd[0]=(dj-cx)/fx;
//   vd[1]=(di-cy)/fy;
//   vd[2]=1.;
//   //printf("vCam= %f\t%f\t%f\n",vd[0],vd[1],vd[2]);
//   Plane2Earth(EA[0],EA[1],EA[2],vd[0],vd[1],vd[2]);//trasformazione di vd nel rif. del drone
//   //printf("vLab= %f\t%f\t%f\n",Vservice[0],Vservice[1],Vservice[2]);
//   t=-Pc[0]/Vservice[0];//form X=Pc[0]+Vservice[0]*t=0
//   XYZ[0]=0.;
//   XYZ[1]=Pc[1]+Vservice[1]*t;
//   XYZ[2]=Pc[2]+Vservice[2]*t;
//   distObs=0.;
//   for(int i=0;i<3;i++)
//       distObs=distObs+(XYZ[i]-Pc[i])*(XYZ[i]-Pc[i]);
//   distObs=sqrt(distObs);
//}



void Plane2Earth(double yaw,double pitch, double roll,double xd, double yd, double zd){
    //From Plane to Earth by rotation matrix Z1Y2X3=Z(-Yaw)Y(-Pitch)X(-Roll)
    // yaw=1; pitch=2; roll=3
    double c1,s1,c2,s2,c3,s3;
    c1=cos(-yaw);
    s1=sin(-yaw);
    c2=cos(-pitch);
    s2=sin(-pitch);
    c3=cos(-roll);
    s3=sin(-roll);
    Vservice[0]=c1*c2*xd +(c1*s2*s3-c3*s1)*yd +(s1*s3+c1*c3*s2)*zd;
    Vservice[1]=c2*s1*xd +(c1*c3+s1*s2*s3)*yd +(c3*s1*s2-c1*s3)*zd;
    Vservice[2]=  -s2*xd +           c2*s3*yd +           c2*c3*zd;
}



void Earth2Plane(double yaw,double pitch, double roll,double xd, double yd, double zd){
    //From Earth to Plane by rotation matrix X1Y2Z3=X(Roll)Y(Pitch)Z(Yaw)
    // yaw=3; pitch=2; roll=1
    double c1,s1,c2,s2,c3,s3;
    c1=cos(roll);
    s1=sin(roll);
    c2=cos(pitch);
    s2=sin(pitch);
    c3=cos(yaw);
    s3=sin(yaw);
    Vservice[0]=           c2*c3*xd            -c2*s3*yd      +s2*zd;
    Vservice[1]=(c1*s3+c3*s1*s2)*xd +(c1*c3-s1*s2*s3)*yd   -c2*s1*zd;
    Vservice[2]=(s1*s3-c1*c3*s2)*xd +(c3*s1+c1*s2*s3)*yd   +c1*c2*zd;
}



void xyz2pxpy(int iCam, double X, double Y, double Z){ //compute pxX pxY given X,Y,Z,Pc[3],EA[3]
    double v[3],vr[3],fx,fy;
    if(iCam==1){
        unitV12(Pc_1[0],Pc_1[1],Pc_1[2],X,Y,Z);//versore target -> drone
        fx=fx_1;
        fy=fy_1;
    }
    else{
        unitV12(Pc_2[0],Pc_2[1],Pc_2[2],X,Y,Z);//versore target -> drone
        fx=fx_2;
        fy=fy_2;
    }
    v[0]=Vservice[0];
    v[1]=Vservice[1];
    v[2]=Vservice[2];
    //printf("vLab= %f\t%f\t%f\n",v[0],v[1],v[2]);
    if(iCam==1)
        Earth2Plane(EA_1[0],EA_1[1],EA_1[2],v[0],v[1],v[2]);//trasformazione di v nel rif. del drone
    else
        Earth2Plane(EA_2[0],EA_2[1],EA_2[2],v[0],v[1],v[2]);//trasformazione di v nel rif. del drone
    //printf("vCam= %f\t%f\t%f\n",Vservice[0],Vservice[1],Vservice[2]);
    vr[0]=Vservice[0];
    vr[1]=Vservice[1];
    vr[2]=Vservice[2];
    pxX=vr[0]*fx/vr[2];
    pxY=vr[1]*fy/vr[2];
}


void PxPyDistorted(int iCam,double xTrue, double yTrue){
    double k1,k2,k3,p1,p2;
    if(iCam==1){
        k1=distCoeffs_1.at<double>(0);
        k2=distCoeffs_1.at<double>(1);
        p1=distCoeffs_1.at<double>(2);
        p2=distCoeffs_1.at<double>(3);
        k3=distCoeffs_1.at<double>(4);
    }
    else{
        k1=distCoeffs_2.at<double>(0);
        k2=distCoeffs_2.at<double>(1);
        p1=distCoeffs_2.at<double>(2);
        p2=distCoeffs_2.at<double>(3);
        k3=distCoeffs_2.at<double>(4);
    }
    double r2=xTrue*xTrue+yTrue*yTrue;
    pxX=xTrue*(1.+k1*r2+k2*pow(r2,2.)+k3*pow(r2,3.))+2.*p1*xTrue*yTrue+p2*(r2+2.*xTrue*xTrue);
    pxY=yTrue*(1.+k1*r2+k2*pow(r2,2.)+k3*pow(r2,3.))+p1*(r2+2.*yTrue*yTrue)+2.*p2*xTrue*yTrue;
}



void unitV12(double x1, double y1, double z1, double x2, double y2, double z2){
    // calculates the unit vector P1->P2
    double sum=0.;
    Vservice[0]=x2-x1;
    Vservice[1]=y2-y1;
    Vservice[2]=z2-z1;
    for(int i=0;i<3;i++)
        sum=sum+Vservice[i]*Vservice[i];
    sum=sqrt(sum);
    for(int i=0;i<3;i++)
        Vservice[i]=Vservice[i]/sum;
}


void sumV12(double x1, double y1, double z1, double x2, double y2, double z2){
    // calculates the unit vector related to the sum of the two vectors V1+V2
    double sum=0.;
    Vservice[0]=x2+x1;
    Vservice[1]=y2+y1;
    Vservice[2]=z2+z1;
    for(int i=0;i<3;i++)
        sum=sum+Vservice[i]*Vservice[i];
    sum=sqrt(sum);
    for(int i=0;i<3;i++)
        Vservice[i]=Vservice[i]/sum;
}


void pxColor(double val){
    //compute the pixel color to build a color map
    int Blue=0;
    int Green=0;
    int Red=0;
    if(val==0.){
        Blue=255;
        Green=0;
        Red=0;
    }
    else if(val>0. && val< 0.25){
        Blue=255;
        Green=int(255.*4.*val+0.5);
        Red=0;
    }
    else if(val>= 0.25 && val< 0.50){
        Blue=int(255.*(1.-4.*(val-0.25)));
        Green=255;
        Red=0;
    }
    else if(val>= 0.50 && val< 0.75){
        Blue=0;
        Green=255;
        Red=int(255.*4.*(val-0.50));
    }
    else if(val>= 0.75 && val<=1.0){
        Blue=0;
        Green=int(255.*(1.-4.*(val-0.75)));
        Red=255;
    }
    else if(val>1){
        Blue=255;
        Green=255;
        Red=255;
    }

    Vservice[0]=Blue;
    Vservice[1]=Green;
    Vservice[2]=Red;
}


int getPosition(){//return x motor current
    int x;
    QString command=part1cmd+"get_position_mm.cmd";
    printf("command= %s\n",command.toStdString().c_str());
    system(command.toStdString().c_str());
    waitKey(100);
    FILE *fd=fopen("/mnt/data/getPosition.txt","r");
    if(fd==NULL)
        return(-1);
    fscanf(fd,"%d",&x);
    fclose(fd);
    return(x);
}


int NINT(double x){
    int n;
    if(x>=0.)
        n=int(x+0.5);
    else
        n=int(x-0.5);
    return(n);
}
