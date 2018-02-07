//
// detect.cpp : Detect integrated circuits in printed circuit board (PCB) images.
//
// Based on skeleton code by D. Crandall, Spring 2018
//
// PUT YOUR NAMES HERE
//
//

#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "fft.h"
#include <vector>
#include "SImage.h"

using namespace std;

struct corn{
    int x;
    int y;
    float value;
};

struct boxes{
    int Ox;
    int Oy;
    int Dx;
    int Dy;
    float value;
};


// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc.
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose.

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);

    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedBox class may be helpful!
//  Feel free to modify.
//
class DetectedBox {
public:
  int row, col, width, height;
  double confidence;
};

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<DetectedBox> &ics)
{
  ofstream ofs(filename.c_str());

  for(vector<DetectedBox>::const_iterator s=ics.begin(); s != ics.end(); ++s)
    ofs << s->row << " " << s->col << " " << s->width << " " << s->height << " " << s->confidence << endl;
}

// Function that outputs a visualization of detected boxes
void  write_detection_image(const string &filename, const vector<DetectedBox> &ics, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];

  for(int p=0; p<3; p++)
    {
      output_planes[p] = input;
      for(vector<DetectedBox>::const_iterator s=ics.begin(); s != ics.end(); ++s)
	overlay_rectangle(output_planes[p], s->row, s->col, s->row+s->height-1, s->col+s->width-1, p==2?255:0, 2);
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}


// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const SDoublePlane &input, SDoublePlane &fft_real, SDoublePlane &fft_imag)
{
    fft_real = input;
    fft_imag = SDoublePlane(input.rows(), input.cols());

    FFT_2D(1, fft_real, fft_imag);
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const SDoublePlane &input_real, const SDoublePlane &input_imag, SDoublePlane &output_real)
{
    output_real = input_real;
    SDoublePlane output_imag = input_imag;

    FFT_2D(0, output_real, output_imag);
}

// The rest of these functions are incomplete. These are just suggestions to
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
    SDoublePlane output(input.rows(), input.cols());

    int offset = row_filter.cols();
    SDoublePlane tempImage(input.rows()+offset,input.cols()+offset);

    int a = offset/2;

    for(int i=0;i<(output.rows());i++){
        for(int j=0;j<(output.cols());j++) {
            tempImage[i+a][j+a] = input[i][j];
        }
    }

    //cols
    for(int i=0;i<(tempImage.rows()-offset);i++){
        for(int j=0;j<(tempImage.cols()-offset);j++){

            int value = 0;
            for(int q=0;q<col_filter.rows();q++){
                value += tempImage[i+q][j]*col_filter[q][0];
            }
            output[i][j]=value;
        }
    }


    for(int i=0;i<(output.rows());i++){
        for(int j=0;j<(output.cols());j++) {
            tempImage[i+a][j+a] = output[i][j];
        }
    }
    //rows
    for(int i=0;i<(tempImage.rows()-offset);i++){
        for(int j=0;j<(tempImage.cols()-offset);j++){

            int value = 0;

            for(int q=0;q<col_filter.cols();q++){
                value += tempImage[i][j+q]*row_filter[0][q];
            }

            output[i][j]=value;
        }
    }
    return output;
}

// Convolve an image with a  convolution kernel
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
    SDoublePlane output(input.rows(), input.cols());

    int offset = filter.cols();
    int shift = offset/2;


    SDoublePlane tempImage(input.rows()+offset,input.cols()+offset);

    for(int i=0;i<(output.rows());i++){
        for(int j=0;j<(output.cols());j++) {
            tempImage[i+shift][j+shift] = input[i][j];
        }
    }
    // Convolution code here
    for(int i=0;i<(tempImage.rows()-offset);i++){
        for(int j=0;j<(tempImage.cols()-offset);j++){

            int value = 0;

            for(int p=0;p<filter.rows();p++){
                for(int q=0;q<filter.cols();q++){
                    value += tempImage[i+p][j+q]*filter[p][q];
                }
            }
            output[i][j]=value;
        }
    }
  return output;
}


// Apply a sobel operator to an image, returns the result
SDoublePlane sobel_gradient_filter(const SDoublePlane &input,int thresh=0)
{
    SDoublePlane output(input.rows(), input.cols());

    //Sobel Filters
    SDoublePlane Xfilter(3,3);

    Xfilter[0][0] = 1;
    Xfilter[0][1] = 2;
    Xfilter[0][2] = 1;

    Xfilter[2][0] = -1;
    Xfilter[2][1] = -2;
    Xfilter[2][2] = -1;

    SDoublePlane Gx = convolve_general(input, Xfilter);

    SDoublePlane Yfilter(3,3);

    Yfilter[0][0] = 1;
    Yfilter[1][0] = 2;
    Yfilter[2][0] = 1;

    Yfilter[0][2] = -1;
    Yfilter[1][2] = -2;
    Yfilter[2][2] = -1;

    SDoublePlane Gy = convolve_general(output, Yfilter);

    int max = 0;

    //calculate magnitude
    for(int i=0;i<input.rows();i++){
        for(int j=0;j<input.cols();j++){
            int mag = sqrt(pow(Gx[i][j],2)+pow(Gy[i][j],2));
            if(max < mag) {
                max = mag;
            }
            output[i][j]=mag;
        }
    }

    //Suppress with threshold
    for(int i=0;i<input.rows();i++){
        for(int j=0;j<input.cols();j++){
            if(output[i][j] < thresh) {
                output[i][j] = 0;
            }
        }
    }
    return output;
}

//Fuction to draw the lines
SDoublePlane drawlines(const SDoublePlane &input,int r,double th){

    for(int i=0;i<input.rows();i++) {
        for (int j = 0; j < input.cols(); j++) {

            if(r == int(i*cos(th)+j*sin(th))){
                input[i][j] = 200;

            }
        }
    }

    return input;

}

//Hough transform for line detection
SDoublePlane hough(const SDoublePlane &input, const SDoublePlane &orig,int numlines){

    int rows = input.rows();
    int cols = input.cols();
    int A = pow(rows,2) + pow(cols,2);

    int Pmax = abs(sqrt(A));

    int Pmin = 1;
    float pi = 3.14159;

    int thetaMin = 0;
    int thetaMax = 360;

    SDoublePlane output(Pmax*2, thetaMax);

    SDoublePlane final(rows, cols);

    for(int i=0;i<input.rows();i++){
        for(int j=0;j<input.cols();j++){

            if(input[i][j] > 0){
                for(int t=0;t<=thetaMax;t++){
                    double r = 0;
                    if(t!=0){
                        r = t*pi/180;
                    }else{
                        r = 0;
                    }

                    int P = int((i*cos(r))+(j*sin(r)));
                    output[Pmax+P][t] +=1;
                }
            }
        }
    }

    int lines0[100];
    int line0=0;
    int lines90[100];
    int line90=0;
    for(int z=0;z<numlines;z++) {

        int max = 0;
        int indi = 0;
        int indj = 0;


        for (int i = 0; i < output.rows(); i++) {
            for (int j = 0; j < output.cols(); j++) {
                if (output[i][j] > max) {
                    max = output[i][j];
                    indi = i;
                    indj = j;
                }
            }
        }
        int r = indi - Pmax;
        double th = indj * pi / 180;

        for (int i = 0; i < 50; i++) {
            for (int j = 0; j < 50; j++) {
                output[indi + i][indj + j] = 0;
                output[indi - i][indj - j] = 0;
            }
        }

        //normalize vertical lines
        if ((th > 1.54 & th < 1.6) || (th > 4.69 & th < 4.75)) {
            r = abs(r);
            if (line90 == 0) {
                lines90[0] = r;
                line90++;
            } else {
                int flag = 0;
                for (int q = 0; q < line90; q++) {
                    if (lines90[q] < r + 10 && lines90[q] > r - 10) {
                        lines90[q] = (lines90[q] + r) / 2;
                        flag = 1;
                        break;
                    }
                }
                if (flag == 0) {
                    lines90[line90] = r;
                    line90++;
                }
            }
        }

        //normalise horizontal lines
        if ((th > 3.1 & th < 3.19) || (th > -0.05 & th < 0.05)) {
            r = abs(r);
            if (line0 == 0) {
                lines0[0] = r;
                line0++;
            } else {
                int flag = 0;
                for (int q = 0; q < line0; q++) {
                    if (lines0[q] < r + 10 && lines0[q] > r - 10) {
                        lines0[q] = (lines0[q] + r) / 2;
                        flag = 1;
                        break;
                    }
                }
                if (flag == 0) {
                    lines0[line0] = r;
                    line0++;
                }
            }
        }
    }

    //draw vertical
    for(int i=0;i<line90;i++){
        final=drawlines(final,lines90[i],1.5708);
    }
    //draw horizontal
    for(int i=0;i<line0;i++){
        final=drawlines(final,lines0[i],0);
    }

    return final;
}

//Harris Corner Detector
SDoublePlane corners(const SDoublePlane &input,int W,const SDoublePlane &inputImage){
    SDoublePlane output(input.rows(), input.cols());

    SDoublePlane final(input.rows(), input.cols());


    SDoublePlane Ix(input.rows(), input.cols());
    SDoublePlane Iy(input.rows(), input.cols());

    for(int i=0;i<input.rows()-1;i++){
        for(int j=0;j<input.cols()-1;j++){
            Ix[i][j] = input[i][j] - input[i][j+1];
            Iy[i][j] = input[i][j] - input[i+1][j];
        }
    }
    for(int i=0;i<input.rows()-W;i++){
        for(int j=0;j<input.cols()-W;j++){

            int A=0;
            int B=0;
            int C=0;
            for(int p=0;p<W;p++){
                for(int q=0;q<W;q++){
                    A += Ix[i+p][j+q] * Ix[i+p][j+q];
                    B += Ix[i+p][j+q] * Iy[i+p][j+q];
                    C += Iy[i+p][j+q] * Iy[i+p][j+q];
                }
            }
            if( A==C   && A!=0 && C!=0){
                output[i][j] = 200;
            }
        }
    }

    return output;
}

//Finds the factor to square the image to
int findfactor(int dim){
    if(dim < 32){
        return 32;
    }else if(dim < 64){
        return 64;
    }else if(dim < 128){
        return 128;
    }else if(dim < 256){
        return 256;
    }else if(dim < 512){
        return 512;
    }

    return 0;
}

//Pads to make square matrix
SDoublePlane padd(const SDoublePlane &input, int size){

    int Wr = (size - input.rows())/2;
    int Wc = (size - input.cols())/2;

    SDoublePlane output(size,size);
    for(int i=0;i<input.rows();i++){
        for(int j=0;j<input.cols();j++){
            output[i+Wr][j+Wc] = input[i][j];
        }
    }

    return output;
}

//chooses ground truth image based on size of test image
SDoublePlane getGTImage(int dim){
    if(dim == 32){
        return SImageIO::read_png_file("test_32.png");
    }
    if(dim == 128){
        return SImageIO::read_png_file("test_128.png");
    }
    if(dim == 256){
        return SImageIO::read_png_file("test_256.png");
    }
    if(dim == 512){
        return SImageIO::read_png_file("test_512.png");
    }
    if(dim == 64){
        return SImageIO::read_png_file("test_64.png");
    }

    return SImageIO::read_png_file("test_512.png");
}


//Squares input image or ground truth image and finds fft.
SDoublePlane makesquare(corn orig,corn next,SDoublePlane &sub_image,int gt=0){

    int X = next.x - orig.x;
    int Y = next.y - orig.y;

    int Fx = findfactor(X);
    int Fy = findfactor(Y);

    if(Fx==0 || Fy==0){
        return sub_image;
    }

    SDoublePlane I(X,Y);

    if(gt == 0) {
        for (int i = orig.x; i < X + orig.x; i++) {
            for (int j = orig.y; j < Y + orig.y; j++) {
                I[i - orig.x][j - orig.y] = sub_image[i][j];
            }
        }
    }

    //square the image and fft
    if(Fx >= Fy) {
        SDoublePlane out(Fx,Fx);

        if(gt==1){
            I = getGTImage(Fx);
        }

        out = padd(I, Fx);

        SDoublePlane real(Fx,Fx);
        SDoublePlane imag(Fx,Fx);
        fft(out,real,imag);

        return real;
    }
    else{
        SDoublePlane out(Fy,Fy);

        if(gt==1){
            I = getGTImage(Fy);
        }

        out = padd(I,Fy);

        SDoublePlane real(Fx,Fx);
        SDoublePlane imag(Fx,Fx);
        fft(out,real,imag);

        return real;
    }
}


//Cross reference Two images to create a score.
float crossref(SDoublePlane &source,SDoublePlane &GT_real){

    int r = source.rows();

    float value = 0;
    float overall = 0;

    for(int i=0;i<r;i++){
        for(int j=0;j<r;j++){
            if(source[i][j]>0 && GT_real[i][j]>0) {
                value += source[i][j];
            }
            overall += GT_real[i][j];
        }
    }
    return value/overall;
}

// This main file just outputs a few test images. You'll want to change it to do
//  something more interesting!
int main(int argc, char *argv[])
{
    if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }

    string input_filename(argv[1]);
    SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());

    // test step 2 by applying mean filters to the input image
    SDoublePlane filter(3,3);
    for(int i=0;i<filter.rows();i++){
        for(int j=0;j<filter.cols();j++){
            filter[i][j] = 0.1111;
        }
    }

    //Basic pipeline
    //edge detection
    SDoublePlane sobel_edges = sobel_gradient_filter(input_image,400);

    //hough transform
    SDoublePlane output_image = hough(sobel_edges,input_image,200);

    // harris corner
    output_image = corners(output_image,1,sobel_edges);

    //find list of corners
    corn cor[10000];
    int c=0;
    for(int i=0;i<output_image.rows();i++){
        for(int j=0;j<output_image.cols();j++){
            if(output_image[i][j] > 199){
                cor[c].x= i;
                cor[c].y =j;
                c++;
            }
        }
    }

    boxes B[1000];
    int b=0;
    float Vmax = 0;

    for(int i=0;i<c;i++){
        corn origin;
        origin.x = cor[i].x;
        origin.y = cor[i].y;

        float value = 0;
        corn max;
        max.x = 1000;
        max.y = 0;

        for(int j=0;j<c;j++){
            if(origin.x < cor[j].x && origin.y < cor[j].y){

                SDoublePlane img = makesquare(origin,cor[j],input_image,0);

                if(img.rows() != img.cols()){
                    continue;
                }

                SDoublePlane GT = makesquare(origin,cor[j],input_image,1);

                float v = crossref(img,GT);

                if(v > value-0.5){
                    int ox = origin.x-max.x;
                    int oy = origin.y-max.y;
                    int vx = origin.x-cor[j].x;
                    int vy = origin.y-cor[j].y;

                    if(v < value+0.5){
                        if((ox*oy) > vx*vy){
                            max.x = cor[j].x;
                            max.y = cor[j].y;
                            value = v;
                        }
                    }else{
                        max.x = cor[j].x;
                        max.y = cor[j].y;
                        value = v;
                    }

                }
            }
        }
        if(int(value) != 0) {
            B[b].Ox = origin.x;
            B[b].Oy = origin.y;
            B[b].Dx = max.x - origin.x;
            B[b].Dy = max.y - origin.y;
            B[b].value = value;
            b++;

            if(Vmax < value){
                Vmax = value;
            }
        }
    }


    DetectedBox s;
    vector<DetectedBox> ics;

    //Write boxes
    vector<DetectedBox> dummy;

    //score threshold:
    float t=0.5;

    for (int i = 0; i < b; i++) {
        if (B[i].value > (Vmax-t)) {
            s.row = B[i].Ox;
            s.col = B[i].Oy;
            s.width = B[i].Dy;
            s.height = B[i].Dx;
            s.confidence = B[i].value/Vmax;
            ics.push_back(s);
        }
    }


    //Answer to question 3
    write_detection_image("edges.png", dummy, sobel_edges);

    //Answer to detection
    write_detection_txt("detected.txt", ics);
    write_detection_image("detected.png", ics, input_image);

}
