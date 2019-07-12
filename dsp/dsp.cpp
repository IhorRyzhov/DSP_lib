#include "dsp.h"


void DftReal(const double *timeSamples, const uint32_t timeSamplesSize, double *freqRE, double *freqIMG)
{
	for (uint32_t k = 0; k < timeSamplesSize / 2 + 1; k++)
	{
		freqRE[k] = 0;
		freqIMG[k] = 0;

		for (uint32_t i = 0; i < timeSamplesSize; i++)
		{
			double arg = (_2Pi * k * i) / timeSamplesSize;

			freqRE[k] += timeSamples[i] * cos(arg);
			freqIMG[k] += timeSamples[i] * sin(arg);
		}

		freqIMG[k] = -freqIMG[k];
		freqIMG[k] *= 2;
		freqRE[k] *= 2;
		freqIMG[k] /= timeSamplesSize;
		freqRE[k] /= timeSamplesSize;
	}
}

void iDftReal(const double *freqSamplesRE, const double *freqSamplesIMG, const LONG freqSamplesSize, double *time)
{
	double *_freqSamplesRE = new double[freqSamplesSize];
	double *_freqSamplesIMG = new double[freqSamplesSize];

	for (LONG i = 0; i < freqSamplesSize; i++)
	{
		_freqSamplesRE[i] = freqSamplesRE[i];
		_freqSamplesIMG[i] = freqSamplesIMG[i];
	}

	LONG N = (freqSamplesSize - 1) * 2;

	_freqSamplesRE[0] /= 2;
	_freqSamplesRE[freqSamplesSize - 1] /= 2;

	for (LONG i = 0; i < N; i++)
	{
		time[i] = 0;

		for (LONG k = 0; k < freqSamplesSize; k++)
		{
			double arg = (2.0 * Pi * k * i) / N;

			time[i] += _freqSamplesRE[k] * cos(arg) - _freqSamplesIMG[k] * sin(arg);
		}
	}

	delete [] _freqSamplesRE;
	delete [] _freqSamplesIMG;
}

void RectToPolar(const float *RE, const float *IMG, const uint32_t size, float *result)
{
	for (uint32_t i = 0; i < size; i++)
	{
		result[i] = sqrt(RE[i] * RE[i] + IMG[i] * IMG[i]);
	}
}

void Normalise(float *RE, const uint32_t size)
{
	for (uint32_t i = 0; i < size; i++)
	{
		RE[i] /= size;
	}
}

void DftComplex(const double *timeSamplesRE, LONG timeSamplesSize, double *freqRE, double *freqIMG)
{
	for (LONG k = 0; k < timeSamplesSize; k++)
	{
		freqRE[k] = 0;
		freqIMG[k] = 0;

		for (LONG n = 0; n < timeSamplesSize; n++)
		{
			double arg = (2.0 * Pi * k * n) / timeSamplesSize;
			freqRE[k] += timeSamplesRE[n] * cos(arg);
			freqIMG[k] += timeSamplesRE[n] * sin(arg);
		}

		freqIMG[k] /= timeSamplesSize;
		freqRE[k] /= timeSamplesSize;
	}
}

void iDftComplex(const double *freqSamplesRE, const double *freqSamplesIMG, const LONG freqSamplesSize, double *timeRE)
{
	double _timeRE = 0;
	double _timeIMG = 0;

	for (LONG k = 0; k < freqSamplesSize; k++)
	{
		for (LONG n = 0; n < freqSamplesSize; n++)
		{
			double arg = (2.0 * Pi * k * n) / freqSamplesSize;

			_timeRE += (freqSamplesRE[n] * cos(arg) + freqSamplesRE[n] * sin(arg));
			_timeIMG += (freqSamplesIMG[n] * cos(arg) - freqSamplesIMG[n] * sin(arg));
		}	
		timeRE[k] = _timeRE - _timeIMG;
		_timeRE = 0;
		_timeIMG = 0;
	}
}

void Convolution(const double* x, const LONG xSize, const double *h, const LONG hSize, double *result)
{
	for (long i = 0; i < xSize + hSize - 1; i++)
	{
		result[i] = 0;

		for (int j = 0; j < hSize; j++)
		{
			if ((i - j < xSize) && (i - j >= 0))
			{
				result[i] += h[j] * x[i - j];
			}
		}
	}
}

LONG bitReverse(const LONG index, const LONG sizeOfBits)
{
	LONG _index = index;
	LONG mask = 0;
	LONG size = sizeOfBits >> 1;

	for (LONG i = 0; i < size; ++i)
	{
		mask = sizeOfBits - i - 1;
		if (_index&(1 << i) && (!(_index&(1 << mask))) || _index&(1 << mask) && (!(_index&(1 << i))))
		{
			_index ^= ((1 << i) | (1 << mask));
		}
	}
	return _index;
}

void fft(const double *timeSamplesRE, const LONG timeSamplesSize, LONG log, double *freqRE, double *freqIMG)
{
	double temp1RE, temp1IMG, temp2RE, temp2IMG, temp3RE, temp3IMG;
	LONG index;
	LONG numberBf = timeSamplesSize / 2; // количество бабочек на определенном этапе. Для 16 точек - 8 двойных бабочки
	LONG sizeBf = timeSamplesSize / numberBf - 1; // ширина бабочек. Начальное значение - 1
	LONG indexL, indexH;
	double arg;

	for (LONG stepFft = 0; stepFft < log; ++stepFft) // 0 ... 4
	{
		index = 0;

		for (LONG i = 0; i < numberBf; ++i) // 8 .. 4 .. 2 .. 1
		{
			for (LONG j = 0; j < sizeBf; ++j) // 1 .. 2 .. 4 .. 8
			{
				indexL = index + j;
				indexH = indexL + sizeBf;

				if (stepFft == 0) // на нулевом шаге производим 2-точечные бпф с инвертированым порядком бит
				{
					temp1RE = timeSamplesRE[bitReverse(indexL, log)];
					temp2RE = timeSamplesRE[bitReverse(indexH, log)];
					freqRE[indexL] = temp1RE + temp2RE;
					freqRE[indexH] = temp1RE - temp2RE;
					freqIMG[indexL] = 0;
					freqIMG[indexH] = 0;
				}
				else
				{
					temp1RE = freqRE[indexL];
					temp1IMG = freqIMG[indexL];					
					temp2RE = freqRE[indexH];
					temp2IMG = freqIMG[indexH];

					arg = _2Pi * j / (sizeBf << 1);

					temp3RE = temp2RE * cos(arg) - temp2IMG * sin(arg);
					temp3IMG = temp2IMG * cos(arg) + temp2RE * sin(arg);

					freqRE[indexL] = temp1RE + temp3RE;
					freqIMG[indexL] = temp1IMG + temp3IMG;

					freqRE[indexH] = temp1RE - temp3RE;
					freqIMG[indexH] = temp1IMG - temp3IMG;

				}
			}
			index += (sizeBf << 1);
		}		
		sizeBf = sizeBf << 1; // sizeBf *= 2;
		numberBf = numberBf >> 1; // numberBf /= 2;
	}
}

void ifft(const double *freqSamplesRE, const double *freqSamplesIMG, const LONG freqSamplesSize, LONG log, double *timeRE, double *timeIMG)
{
	double temp1RE, temp1IMG, temp2RE, temp2IMG, temp3RE, temp3IMG;
	LONG index;
	LONG numberBf = freqSamplesSize / 2; // количество бабочек на определенном этапе. Для 16 точек - 8 двойных бабочки
	LONG sizeBf = freqSamplesSize / numberBf - 1; // ширина бабочек. Начальное значение - 1
	LONG indexL, indexH;
	double arg;

	for (LONG stepFft = 0; stepFft < log; ++stepFft) // 0 ... 4
	{
		index = 0;

		for (LONG i = 0; i < numberBf; ++i) // 8 .. 4 .. 2 .. 1
		{
			for (LONG j = 0; j < sizeBf; ++j) // 1 .. 2 .. 4 .. 8
			{
				indexL = index + j;
				indexH = indexL + sizeBf;

				if (stepFft == 0) // на нулевом шаге производим 2-точечные бпф с инвертированым порядком бит
				{
					temp1RE = freqSamplesRE[bitReverse(indexL, log)];
					temp1IMG = freqSamplesIMG[bitReverse(indexL, log)];
					temp2RE = freqSamplesRE[bitReverse(indexH, log)];
					temp2IMG = freqSamplesIMG[bitReverse(indexH, log)];

					timeRE[indexL] = temp1RE + temp2RE;
					timeIMG[indexL] = temp1IMG + temp2IMG;
					timeRE[indexH] = temp1RE - temp2RE;
					timeIMG[indexH] = temp1IMG - temp2IMG;
				}
				else
				{
					temp1RE = timeRE[indexL];
					temp1IMG = timeIMG[indexL];
					temp2RE = timeRE[indexH];
					temp2IMG = timeIMG[indexH];

					arg = -_2Pi * j / (sizeBf << 1);

					temp3RE = temp2RE * cos(arg) - temp2IMG * sin(arg);
					temp3IMG = temp2IMG * cos(arg) + temp2RE * sin(arg);

					timeRE[indexL] = temp1RE + temp3RE;
					timeIMG[indexL] = temp1IMG + temp3IMG;

					timeRE[indexH] = temp1RE - temp3RE;
					timeIMG[indexH] = temp1IMG - temp3IMG;

				}
			}
			index += (sizeBf << 1);
		}
		sizeBf = sizeBf << 1; // sizeBf *= 2;
		numberBf = numberBf >> 1; // numberBf /= 2;
	}
}

void fastConvolution(const double* x, const double *h, double *resultRE, double *resultIMG, const LONG resultConvolutionSize, const LONG log)
{
	double *_hRE = new double[resultConvolutionSize];
	double *_hIMG = new double[resultConvolutionSize];
	double RE, IMG;

	fft(x, resultConvolutionSize, log, resultRE, resultIMG);
	fft(h, resultConvolutionSize, log, _hRE, _hIMG);

	for (LONG i = 0; i < resultConvolutionSize; i++)
	{
		RE = resultRE[i] * _hRE[i] - resultIMG[i] * _hIMG[i];
		IMG = resultIMG[i] * _hRE[i] + resultRE[i] * _hIMG[i];
		
		_hRE[i] = RE;
		_hIMG[i] = IMG;
	}

	ifft(_hRE, _hIMG, resultConvolutionSize, log, resultRE, resultIMG);
	for (LONG i = 0; i < resultConvolutionSize; i++)
	{
		resultRE[i] /= resultConvolutionSize;
		resultIMG[i] /= resultConvolutionSize;
	}

	delete[] _hRE;
	delete[] _hIMG;
}


