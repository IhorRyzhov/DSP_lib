#pragma once


#include <windows.h>
#include <cmath> 
#include <stdint.h>

#define _2Pi 6.28318
#define Pi 3.14159

/* ������������ ���������� �������������� �����

*timeSamples     - ������ ������� �� ��������� �������;
timeSamplesSize - ������ ������� ������� �� ��������� �������;
*freqRE          - ������ �������������� ��������, ������ ������� = size / 2 + 1;
*freqIMG         - ������ ������ ��������, ������ ������� = size / 2 + 1;

*/
void DftReal(const double *timeSamples, const uint32_t timeSamplesSize, double *freqRE, double *freqIMG);


/* ������������ �������� ���������� �������������� �����

*freqSamplesRE   - ������ �������������� ��������;
*freqSamplesIMG  - ������ ������ ��������;
freqSamplesSize  - ������ ������� ������� � ��������� �������;
*time            - ������ ������� �� ��������� �������;

*/
void iDftReal(const double *freqSamplesRE, const double *freqSamplesIMG, const LONG freqSamplesSize, double *time);


/* ������� �� ������������� � �������� ������� ���������

RE     - ������ ��������������� ��������;
IMG    - ������ ������ ��������;
size   - ������ �������� RE, IMG;
result - ������ ������� �������� � �������� �����������;

*/
void RectToPolar(const float *RE, const float *IMG, const uint32_t size, float *result);

void Normalise(float *RE, const uint32_t size);

/* ����������� ���������� �������������� �����

*timeSamplesRE   - ������ �������������� ��������;
timeSamplesSize  - ������ ������� ������� �� ��������� �������;
*freqSamplesRE   - ������ �������������� ��������, ������ ������� = timeSamplesSize;
*freqSamplesIMG  - ������ ������ ��������, ������ ������� = timeSamplesSize;

*/
void DftComplex(const double *timeSamplesRE, const LONG timeSamplesSize, double *freqSamplesRE, double *freqSamplesIMG);


/* ����������� �������� ���������� �������������� �����

*freqSamplesRE   - ������ ������� � ��������� �������;
*freqSamplesIMG  - ������ ������ ��������;
freqSamplesSize - ������ ������� ������� � ��������� �������;
*time            - ������ ������� �� ��������� �������;

*/
void iDftComplex(const double *freqSamplesRE, const double *freqSamplesIMG, const LONG freqSamplesSize, double *timeRE);

/* �������� �������

*/
void Convolution(const double* x, const LONG xSize, const double *h, const LONG hSize, double *result);

/* ����������� ���

*/
void fft(const double *timeSamplesRE, const LONG timeSamplesSize, LONG log, double *freqRE, double *freqIMG);

/* ����������� ����

*/
void ifft(const double *freqSamplesRE, const double *freqSamplesIMG, const LONG freqSamplesSize, LONG log, double *timeRE, double *timeIMG);

/* ������� �������� �������

resultConvolutionSize = size(x) + size(h) - 1 � ������ ���� ����� : 128, 256, 512, 1024, ��������: size(x) = 150, size(h) = 107, 150 + 107 - 1 = 256

x, h - ������ ����� ������ resultConvolutionSize : 128, 256, 512, 1024 � ��� ������ ���� ��������� ������ (��� ����� ��� ���)
log - ������� ������ resultConvolutionSize : 7, 8, 9, 10

*/
void fastConvolution(const double* x, const double *h, double *resultRE, double *resultIMG, const LONG resultConvolutionSize, const LONG log);
