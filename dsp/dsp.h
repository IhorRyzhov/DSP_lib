#pragma once


#include <windows.h>
#include <cmath> 
#include <stdint.h>

#define _2Pi 6.28318
#define Pi 3.14159

/* ¬ещественное дискретное преобразование ‘урье

*timeSamples     - массив выборок во временной области;
timeSamplesSize - размер массива выборок во временной области;
*freqRE          - массив действительных значений, размер массива = size / 2 + 1;
*freqIMG         - массив мнимых значений, размер массива = size / 2 + 1;

*/
void DftReal(const double *timeSamples, const uint32_t timeSamplesSize, double *freqRE, double *freqIMG);


/* ¬ещественное обратное дискретное преобразование ‘урье

*freqSamplesRE   - массив действительных значений;
*freqSamplesIMG  - массив мнимых значений;
freqSamplesSize  - размер массива выборок в частотной области;
*time            - массив выборок во временной области;

*/
void iDftReal(const double *freqSamplesRE, const double *freqSamplesIMG, const LONG freqSamplesSize, double *time);


/* ѕереход от пр€моугольной к пол€рной системе координат

RE     - массив действичтельных значений;
IMG    - массив мнимых значений;
size   - размер массивов RE, IMG;
result - массив модулей значений в пол€рных координатах;

*/
void RectToPolar(const float *RE, const float *IMG, const uint32_t size, float *result);

void Normalise(float *RE, const uint32_t size);

/*  омплексное дискретное преобразование ‘урье

*timeSamplesRE   - массив действительных значений;
timeSamplesSize  - размер массива выборок во временной области;
*freqSamplesRE   - массив действительных значений, размер массива = timeSamplesSize;
*freqSamplesIMG  - массив мнимых значений, размер массива = timeSamplesSize;

*/
void DftComplex(const double *timeSamplesRE, const LONG timeSamplesSize, double *freqSamplesRE, double *freqSamplesIMG);


/*  омплексное обратное дискретное преобразование ‘урье

*freqSamplesRE   - массив выборок в частотной области;
*freqSamplesIMG  - массив мнимых значений;
freqSamplesSize - размер массива выборок в частотной области;
*time            - массив выборок во временной области;

*/
void iDftComplex(const double *freqSamplesRE, const double *freqSamplesIMG, const LONG freqSamplesSize, double *timeRE);

/* Ћинейна€ свертка

*/
void Convolution(const double* x, const LONG xSize, const double *h, const LONG hSize, double *result);

/*  омплексное Ѕѕ‘

*/
void fft(const double *timeSamplesRE, const LONG timeSamplesSize, LONG log, double *freqRE, double *freqIMG);

/*  омплексное ќЅѕ‘

*/
void ifft(const double *freqSamplesRE, const double *freqSamplesIMG, const LONG freqSamplesSize, LONG log, double *timeRE, double *timeIMG);

/* Ѕыстра€ линейна€ свертка

resultConvolutionSize = size(x) + size(h) - 1 и должен быть равен : 128, 256, 512, 1024, например: size(x) = 150, size(h) = 107, 150 + 107 - 1 = 256

x, h - должны иметь размер resultConvolutionSize : 128, 256, 512, 1024 и уже должны быть дополнены нул€ми (так нужно дл€ бпф)
log - степень двойки resultConvolutionSize : 7, 8, 9, 10

*/
void fastConvolution(const double* x, const double *h, double *resultRE, double *resultIMG, const LONG resultConvolutionSize, const LONG log);
