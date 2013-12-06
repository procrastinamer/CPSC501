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

struct WAV_HEADER
{
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

int main(int argc, char* argv[])
{
	if (argc == 4)
	{
		cout << "Input: " << argv[1] << "\n"
			<< "Impulse: " << argv[2] << "\n"
			<< "Output: " << argv[3] << "\n";

		cout << "Size of WAV_HEADER: " << sizeof(WAV_HEADER) << "\n\n";

		struct WAV_HEADER header;
		struct WAV_HEADER irheader;
		
		FILE * input_file;
		FILE * impulse_file;
		FILE * output_file;
		
		input_file = fopen(argv[1], "rb");
		impulse_file = fopen(argv[2], "rb");
		output_file = fopen(argv[3], "wb");
		
		if (input_file == NULL || impulse_file == NULL || output_file == NULL)
		{
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

		if (header.num_channels > 2)
		{
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
			
		if (irheader.num_channels > 2)
		{
			cout << "Unsupported channel number: " << irheader.num_channels << "\n";
			return -1;
		}
		
		// Calculate buffer sizes
		size_t bytes_per_sample = header.bits_per_sample / 8;
		size_t N = header.subchunk2_size / (bytes_per_sample * header.num_channels);
		
		cout << "Bytes per sample: " << bytes_per_sample << "\n";
		cout << "N (Number of samples) = " << N << "\n\n";
		
		bytes_per_sample = irheader.bits_per_sample / 8;
		size_t M = irheader.subchunk2_size / (bytes_per_sample * irheader.num_channels);

		int P = N + M - 1;

		cout << "Bytes per sample: " << bytes_per_sample << "\n";
		cout << "M (Number of samples) = " << M << "\n\n";
		
		// initialize buffers
		double * in_buf1 = new double[N];
		double * in_buf2 = new double[N];
		
		double * ir_buf1 = new double[M];
		double * ir_buf2 = new double[M];
		
		double * out_buf1 = new double[P];
		double * out_buf2 = new double[P];
		
		// Read input data
		
		if (bytes_per_sample == 1)
		{
			int8_t sample;
			for (int i = 0; i < N; i++)
			{
				fread(&sample, bytes_per_sample, 1, input_file);				
				in_buf1[i] = (double)sample;
				if (header.num_channels == 2)
				{
					fread(&sample, bytes_per_sample, 1, input_file);
					in_buf2[i] = (double)sample;
				}
			}
		}
		else if (bytes_per_sample == 2)
		{	
			int16_t sample;
			for (int i = 0; i < N; i++)
			{
				fread(&sample, bytes_per_sample, 1, input_file);
				in_buf1[i] = (double)sample;
				if (header.num_channels == 2)
				{
					fread(&sample, bytes_per_sample, 1, input_file);
					in_buf2[i] = (double)sample;
				}
			}
		}
		else
		{
			cout << "Unsupported sample size " << bytes_per_sample << "\n";
		}
		
		// Read IR data
		if (bytes_per_sample == 1)
		{
			double sample;
			for (int i = 0; i < M; i++)
			{
				fread(&sample, bytes_per_sample, 1, impulse_file);
				ir_buf1[i] = sample;
				
				if (irheader.num_channels == 2)
				{
					fread(&sample, bytes_per_sample, 1, impulse_file);
					ir_buf2[i] = sample;
				}
			}
		}
		else if (bytes_per_sample == 2)
		{
			double sample;
			for (int i = 0; i < M; i++)
			{
				fread(&sample, bytes_per_sample, 1, impulse_file);
				ir_buf1[i] = sample;
				
				if (irheader.num_channels == 2)
				{
					fread(&sample, bytes_per_sample, 1, impulse_file);
					ir_buf2[i] = sample;
				}
			}
		}
		else
		{
			cout << "Unsupported sample size " << bytes_per_sample << "\n";
			return -1;
		}
		
		// Convolve
		cout << "Now convolving..." << "\n";
		
		clock_t timer = clock();
		
		for (int n = 0; n < N; n++)
		{
			for (int m = 0; m < M; m++)
			{			
				out_buf1[n + m] += in_buf1[n] * ir_buf1[m];
				if (header.num_channels == 2) out_buf2[n + m] += ir_buf2[m] * in_buf2[n];
			}
		}

		timer = clock() - timer;
		
		cout << "Convolution took " << timer/CLOCKS_PER_SEC << " seconds \n";

		// Write data
		
		fwrite(&header, sizeof(WAV_HEADER), 1, output_file);
		
		if (bytes_per_sample == 1)
		{
			double sample;
			for (int i = 0; i < N; i++)
			{
				sample = out_buf1[i];
				fwrite(&sample, bytes_per_sample, 1, output_file);
				
				if (header.num_channels == 2)
				{
					sample = out_buf2[i];
					fwrite(&sample, bytes_per_sample, 1, output_file);
				}
			}
		}
		else if (bytes_per_sample == 2)
		{
			double sample;
			for (int i = 0; i < N; i++)
			{
				sample = out_buf1[i];
				fwrite(&sample, bytes_per_sample, 1, output_file);

				if (header.num_channels == 2)
				{
					sample = out_buf2[i];
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
