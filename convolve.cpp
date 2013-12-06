// Douglas Lam
// 10053718
// CPSC 501 A4

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <time.h>
using namespace std;

#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

struct WAV_HEADER {
	// Header 1
	char RIFF[4];
	uint32_t chunk_size;
	char WAVE[4];

	// Header 2
	char fmt[4];
	uint32_t subchunk1_size;
	uint16_t audio_format;
	uint16_t num_channels;
	uint32_t sample_rate;
	uint32_t byte_rate;
	uint16_t block_align;
	uint16_t bits_per_sample;
	uint16_t extra;

	// Data
	char data[4];
	uint32_t subchunk2_size;
} __attribute__((packed));

//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
void four1(double data[], int nn, int isign) {
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}

int main(int argc, char* argv[]) {
	if (argc == 4) {
		cout << "Input: " << argv[1] << "\n"
			<< "Impulse: " << argv[2] << "\n"
			<< "Output: " << argv[3] << "\n";

		cout << "Size of WAV_HEADER: " << sizeof(WAV_HEADER) << "\n\n";

		int i;
		int ii;
		
		struct WAV_HEADER header;
		struct WAV_HEADER irheader;
		
		FILE * input_file;
		FILE * impulse_file;
		FILE * output_file;
		
		input_file = fopen(argv[1], "rb");
		impulse_file = fopen(argv[2], "rb");
		output_file = fopen(argv[3], "wb");
		
		if (input_file == NULL || impulse_file == NULL || output_file == NULL) {
			cout << "Cannot open the file.\n";
			return -1;
		}

		// Read input header
		fread(&header, sizeof(WAV_HEADER), 1, input_file);

		cout << "ChunkID: '" << header.RIFF << "'\n"
			<< "ChunkSize: " << header.chunk_size << "\n"
			<< "Format: '" << header.WAVE << "'\n"
			<< "SubChunk1ID: '" << header.fmt << "'\n"
			<< "SubChunk1Size: " << header.subchunk1_size << "\n"
			<< "AudioFormat: " << header.audio_format << "\n"
			<< "NumChannels: " << header.num_channels << "\n"
			<< "SampleRate: " << header.sample_rate << "\n"
			<< "ByteRate: " << header.byte_rate << "\n"
			<< "BlockAlign: " << header.block_align << "\n"
			<< "BitsPerSample: " << header.bits_per_sample << "\n"
			<< "SubChunk2ID: '" << header.data << "'\n"
			<< "SubChunk2Size: " << header.subchunk2_size << "\n\n";

		if (header.num_channels > 2) {
			cout << "Unsupported channel number: " << header.num_channels << "\n";
		}
				
		// Read IR header
		fread(&irheader, sizeof(WAV_HEADER), 1, impulse_file);

		cout << "ChunkID: '" << irheader.RIFF << "'\n"
			<< "ChunkSize: " << irheader.chunk_size << "\n"
			<< "Format: '" << irheader.WAVE << "'\n"
			<< "SubChunk1ID: '" << irheader.fmt << "'\n"
			<< "SubChunk1Size: " << irheader.subchunk1_size << "\n"
			<< "AudioFormat: " << irheader.audio_format << "\n"
			<< "NumChannels: " << irheader.num_channels << "\n"
			<< "SampleRate: " << irheader.sample_rate << "\n"
			<< "ByteRate: " << irheader.byte_rate << "\n"
			<< "BlockAlign: " << irheader.block_align << "\n"
			<< "BitsPerSample: " << irheader.bits_per_sample << "\n"
			<< "SubChunk2ID: '" << irheader.data << "'\n"
			<< "SubChunk2Size: " << irheader.subchunk2_size << "\n\n";
			
		if (irheader.num_channels > 2) {
			cout << "Unsupported channel number: " << irheader.num_channels << "\n";
			return -1;
		}
		
		// Calculate buffer sizes
		size_t bytes_per_sample = header.bits_per_sample / 8;
		size_t N = header.subchunk2_size / (bytes_per_sample * header.num_channels);
		
		cout << "Bytes per sample: " << bytes_per_sample << "\n";
		cout << "N (Number of samples) = " << N << "\n\n";
		
		size_t ir_bytes_per_sample = irheader.bits_per_sample / 8;
		size_t M = irheader.subchunk2_size / (ir_bytes_per_sample * irheader.num_channels);

		cout << "Bytes per sample: " << ir_bytes_per_sample << "\n";
		cout << "M (Number of samples) = " << M << "\n\n";
		
		// Make the size a power of 2
		int size = pow(2, ceil(log(N)/log(2)));
		cout << "Size: " << size << "\n";
		
		// initialize buffers
		double * in_buf1 = new double[size*2];
		double * in_buf2 = new double[size*2];
		
		double * ir_buf1 = new double[size*2];
		double * ir_buf2 = new double[size*2];
		
		double * out_buf1 = new double[size*2];
		double * out_buf2 = new double[size*2];

		for (i = 0; i < size*2; i++) {
			in_buf1[i] = 0.0;
			in_buf2[i] = 0.0;
			ir_buf1[i] = 0.0;
			ir_buf2[i] = 0.0;
			out_buf1[i] = 0.0;
			out_buf2[i] = 0.0;
		}
		
		// Read input data
		double largest = 0.0;
		if (bytes_per_sample == 1) {
			int8_t sample;
			for (i = 0, ii = 0; i < N; i++, ii += 2) {
				fread(&sample, bytes_per_sample, 1, input_file);
				in_buf1[ii] = (double)sample;
	
				if (abs(in_buf1[ii]) >= largest) { largest = abs(in_buf1[ii]); }

				if (header.num_channels == 2) {
					fread(&sample, bytes_per_sample, 1, input_file);
					in_buf2[ii] = (double)sample;

					if (abs(in_buf2[ii]) >= largest) { largest = abs(in_buf2[ii]); }
				}
			}
		}
		else if (bytes_per_sample == 2) {
			int16_t sample;
			for (i = 0, ii = 0; i < N; i++, ii += 2) {
				fread(&sample, bytes_per_sample, 1, input_file);
				in_buf1[ii] = (double)sample;
	
				if (abs(in_buf1[ii]) >= largest) { largest = abs(in_buf1[ii]); }
	
				if (header.num_channels == 2) {
					fread(&sample, bytes_per_sample, 1, input_file);
					in_buf2[ii] = (double)sample;
	
					if (abs(in_buf2[ii]) >= largest) { largest = abs(in_buf2[ii]); }
				}
			}
		}
		else {
			printf("Unsupported sample size %d\n", bytes_per_sample);
			return -1;
		}

		// Read IR data
		if (ir_bytes_per_sample == 1) {
			int8_t sample;
			for (i = 0, ii = 0; i < M; i++, ii += 2) {
				fread(&sample, ir_bytes_per_sample, 1, impulse_file);
				ir_buf1[ii] = (double)sample;
				
				if (abs(ir_buf1[ii]) >= largest) { largest = abs(ir_buf1[ii]); }
				
				if (irheader.num_channels == 2) {
					fread(&sample, ir_bytes_per_sample, 1, impulse_file);
					ir_buf2[ii] = (double)sample;
					
					if (abs(ir_buf2[ii]) >= largest) { largest = abs(ir_buf2[ii]); }
				}
			}
		}
		else if (ir_bytes_per_sample == 2) {
			int16_t sample;
			for (i = 0, ii = 0; i < M; i++, ii += 2) {
				fread(&sample, ir_bytes_per_sample, 1, impulse_file);
				ir_buf1[ii] = (double)sample;

				if (abs(ir_buf1[ii]) >= largest) { largest = abs(ir_buf1[ii]); }

				if (irheader.num_channels == 2) {
					fread(&sample, ir_bytes_per_sample, 1, impulse_file);
					ir_buf2[ii] = (double)sample;

					if (abs(ir_buf2[ii]) >= largest) { largest = abs(ir_buf2[ii]); }
				}
			}
		}
		else {
			printf("Unsupported sample size %d\n", ir_bytes_per_sample);
			return -1;
		}
		
		cout << "Largest: " << largest << "\n";
		
		// Normalize
		for (i = 0; i < size*2; i++) {
			in_buf1[i] /= largest;
			in_buf2[i] /= largest;
			
			ir_buf1[i] /= largest;
			ir_buf2[i] /= largest;
		}

		// Convolve
		cout << "Now convolving..." << "\n";
		
		clock_t timer = clock();
		
		four1(in_buf1-1, size, 1);		
		if (header.num_channels == 2) { four1(in_buf2-1, size, 1); }
		four1(ir_buf1-1, size, 1);
		if (irheader.num_channels == 2) { four1(ir_buf2-1, size, 1); }

		// Complex multiplication
		for (i = 0, ii = 0; i < size; i++, ii += 2) {
			out_buf1[ii] = in_buf1[ii] * ir_buf1[ii] - in_buf1[ii+1] * ir_buf1[ii+1];
			out_buf1[ii+1] = in_buf1[ii] * ir_buf1[ii+1] + in_buf1[ii+1] * ir_buf1[ii];
			if (header.num_channels == 2) {
				out_buf2[ii] = in_buf2[ii] * ir_buf2[ii] - in_buf2[ii+1] * ir_buf2[ii+1];
				out_buf2[ii+1] = in_buf2[ii] * ir_buf2[ii+1] + in_buf2[ii+1] * ir_buf2[ii];
			}
		}
		
		four1(out_buf1-1, size, -1);
		four1(out_buf2-1, size, -1);

		// Scaling
		for (i = 0, ii = 0; i < size; i++, ii+=2) {
					out_buf1[ii] /= size;
					out_buf1[ii+1] /= size;
					out_buf2[ii] /= size;
					out_buf2[ii+1] /= size;
		}
		
		for (i = 0; i < size*2; i++)  
		{
			out_buf1[i] *= largest;
			out_buf1[i] = round(out_buf1[i]);
			out_buf2[i] *= largest;
			out_buf2[i] = round(out_buf2[i]);
		}
		
		timer = clock() - timer;
		
		cout << "Convolution took " << timer/CLOCKS_PER_SEC << " seconds and " << timer << " ticks \n";

		// Write data
		
		fwrite(&header, sizeof(WAV_HEADER), 1, output_file);
		
		if (bytes_per_sample == 1) {
			int8_t sample;
			for (i = 0, ii = 0; i < N; i++, ii +=2) {
				sample = out_buf1[ii];
				fwrite(&sample, bytes_per_sample, 1, output_file);
				if (header.num_channels == 2) {
					sample = out_buf2[ii];
					fwrite(&sample, bytes_per_sample, 1, output_file);
				}
			}
		}
		else if (bytes_per_sample == 2) {
			int16_t sample;
			for (i = 0, ii = 0; i < N; i++, ii +=2) {
				sample = out_buf1[ii];
				fwrite(&sample, bytes_per_sample, 1, output_file);
				if (header.num_channels == 2) {
					sample = out_buf2[ii];
					fwrite(&sample, bytes_per_sample, 1, output_file);
				}
			}
		}
		else
		{
			cout << "Unsupported sample size" << bytes_per_sample << "\n";
			return -1;
		}

		fclose(input_file);
		fclose(impulse_file);
		fclose(output_file);

		return 0;
	}
	else
	{
		cout << "Missing arguments. Please try again with the format: 'convolve input.wav impulse.wav output.wav' \n";
		return -1;
	}
}
