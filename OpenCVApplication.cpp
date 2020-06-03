// OpenCVApplication.cpp : Defines the entry point for the console application.
//

#pragma once

#include "stdafx.h"
#include "common.h"
#include <queue>
#include <stack>
#include <random>
using namespace std;


double c(int i, int n)
{
	return (i == 0) ? sqrt(1.0 / n) : sqrt(2.0 / n);
}


Mat_<double> DCT(Mat_<uchar> src, int n)
{
	Mat_<double> dst(src.rows, src.cols);

	Mat_<int> src1(src.rows, src.cols);

	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			src1(i, j) = src(i, j) - 128;
		}
	}

	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			double temp = 0;

			for (int x = 0; x < n; x++)
			{
				for (int y = 0; y < n; y++)
				{
					temp += src1(x, y) * cos((double)((2 * x + 1) * i * PI) / (2.0 * n)) * cos((double)((2 * y + 1) * j * PI) / (2.0 * n));
				}
			}

			dst(i, j) = c(i, n) * c(j, n) * temp;
		}
	}

	return dst;
}

Mat_<int> quantize(Mat_<double> src, Mat_<uchar> quantizationMatrix)
{
	Mat_<int> dst(src.rows, src.cols);

	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			dst(i, j) = round((double)src(i, j) / quantizationMatrix(i, j));
		}
	}

	return dst;
}


vector <pair<int, int>> runLengthCoding(Mat_<int> src)
{
	vector <pair<int, int>> dst;
	int countZero = 0;

	//cross the matrix in zig-zag
	for (int i = 0; i < 15; i++) {

		if (i % 2 == 1) {
			// down left
			// odd => x increases, y decreases

			//i > 8 => always start from last column, y = 7, x = i-7
			int x = i < 8 ? 0 : i - 7;
			int y = i < 8 ? i : 7;

			while (x < 8 && y >= 0) {

				if (src(x, y) != 0)
				{
					int elem = src(x, y);
					dst.emplace_back(std::make_pair(countZero, elem));
					countZero = 0;
				}
				else
				{
					countZero++;
				}

				x++;
				y--;
			}

		}
		else {
			// up right
			//even => x decreases, y increases

			//i > 8 => always start from botton, x = 7, y = i-7
			int x = i < 8 ? i : 7;
			int y = i < 8 ? 0 : i - 7;
			while (x >= 0 && y < 8) {

				if (src(x, y) != 0)
				{
					int elem = src(x, y);
					dst.emplace_back(std::make_pair(countZero, elem));
					countZero = 0;
				}
				else
				{
					countZero++;
				}

				x--;
				y++;
			}
		}
	}

	return dst;
}

Mat_<int> dequantize(Mat_<int> src, Mat_<uchar> quantizationMatrix)
{
	Mat_<int> dst(src.rows, src.cols);

	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			dst(i, j) = src(i, j) * quantizationMatrix(i, j);
		}
	}

	return dst;
}

Mat_<uchar> IDCT(Mat_<int> src, int n)
{
	Mat_<double> dst(src.rows, src.cols);

	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			double temp = 0;

			//fdct
			for (int x = 0; x < n; x++)
			{
				for (int y = 0; y < n; y++)
				{
					temp += c(x, n) * c(y, n) * src(x, y) * cos((double)((2 * i + 1) * x * PI) / (2.0 * n)) * cos((double)((2 * j + 1) * y * PI) / (2.0 * n));
				}
			}

			dst(i, j) = temp;
		}
	}

	Mat_<uchar> dst1(dst.rows, dst.cols);
	int temp;
	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			temp = round(dst(i, j)) + 128;
			temp = temp > 255 ? 255 : temp;

			dst1(i, j) = temp < 0 ? 0 : temp;
		}
	}

	return dst1;
}

Mat_<int> decodeRle(int rows, int cols, vector<pair<int, int>> data)
{
	Mat_<int> dst(rows, cols);

	int n = data.size();
	int k = 0;

	if (data.size() > 0)
	{
		pair<int, int> crt_value = data[0];
		int zeroes = 0;
		boolean isEob = false;

		for (int i = 0; i < 15; i++) {

			if (i % 2 == 1) {

				int x = i < 8 ? 0 : i - 7;
				int y = i < 8 ? i : 7;

				while (x < 8 && y >= 0) {

					if (!isEob)
					{
						if (zeroes != 0)
						{
							dst(x, y) = 0;
							zeroes--;
						}
						else
						{
							int elem = crt_value.second;
							dst(x, y) = elem;

							k++;

							if (k == n)
							{
								isEob = true;
							}
							else
							{
								crt_value = data[k];
								zeroes = crt_value.first;
							}
						}
					}
					else
					{
						dst(x, y) = 0;
					}

					x++;
					y--;
				}

			}
			else {

				int x = i < 8 ? i : 7;
				int y = i < 8 ? 0 : i - 7;
				while (x >= 0 && y < 8) {

					if (!isEob)
					{
						if (zeroes != 0)
						{
							dst(x, y) = 0;
							zeroes--;
						}
						else
						{
							int elem = crt_value.second;
							dst(x, y) = elem;

							k++;

							if (k == n)
							{
								isEob = true;
							}
							else
							{
								crt_value = data[k];
								zeroes = crt_value.first;
							}
						}


					}
					else
					{
						dst(x, y) = 0;
					}

					x--;
					y++;
				}
			}
		}

		return dst;
	}
	else
	{
		return Mat_<int>(rows, cols, 0);
	}
}

bool isInside(Mat img, int i, int j)
{
	return i >= 0 && i < img.rows && j >= 0 && j < img.cols;
}

vector<vector<pair<int, int>>> compress(Mat_<uchar> src, Mat_<uchar> quantization)
{
	vector<vector<pair<int, int>>> res;
	Mat_<uchar> block(8, 8);

	for (int i = 0; i < src.rows; i = i + 8)
	{
		for (int j = 0; j < src.cols; j = j + 8)
		{
			//get next block starting from  (i, j)
			for (int u = 0; u < 8; u++)
			{
				for (int v = 0; v < 8; v++)
				{
					int i1 = i + u;
					int j1 = j + v;

					if (isInside(src, i1, j1))
					{
						block(u, v) = src(i1, j1);
					}
					else
					{
						block(u, v) = 0;
					}
				}
			}

			//apply DCT
			Mat_<double> afterDCT = DCT(block, 8);

			//quantize
			Mat_<int> afterQuantization = quantize(afterDCT, quantization);

			vector<pair<int, int>> afterRle = runLengthCoding(afterQuantization);

			//rle
			res.push_back(afterRle);
		}
	}


	return res;
}

Mat_<uchar> decompress(int rows, int cols, vector<vector<pair<int, int>>> data, Mat_<uchar> quantization)
{
	Mat_<uchar> dst(rows, cols);

	vector<pair<int, int>> crt_data;
	int k = 0;

	//cross the image from 8x8
	for (int i = 0; i < rows; i = i + 8)
	{
		for (int j = 0; j < cols; j = j + 8)
		{
			if (k < data.size())
			{
				crt_data = data.at(k);

				Mat_<int> decodedImage = decodeRle(8, 8, crt_data);

				//dequantize
				Mat_<int> dequantizedImage = dequantize(decodedImage, quantization);

				//idct
				Mat_<uchar> idctImage = IDCT(dequantizedImage, 8);

				for (int u = 0; u < 8; u++)
				{
					for (int v = 0; v < 8; v++)
					{
						int i1 = i + u;
						int j1 = j + v;

						if (isInside(dst, i1, j1))
						{
							dst(i1, j1) = idctImage(u, v);
						}
					}
				}

				k++;
			}

		}
	}

	return dst;
}

vector<Mat_<uchar>> split(Mat_<Vec3b> src)
{
	vector<Mat_<uchar>> dst;
	Mat_<uchar> Y(src.rows, src.cols);
	Mat_<uchar> U(src.rows, src.cols);
	Mat_<uchar> V(src.rows, src.cols);

	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			Y(i, j) = src(i, j)[0];
			U(i, j) = src(i, j)[1];
			V(i, j) = src(i, j)[2];
		}
	}

	dst.push_back(Y);
	dst.push_back(U);
	dst.push_back(V);

	return dst;
}

Mat_<Vec3b> reassemble(Mat_<uchar> Y, Mat_<uchar> U, Mat_<uchar> V)
{
	Mat_<Vec3b> dst(Y.rows, Y.cols);

	for (int i = 0; i < dst.rows; i++)
	{
		for (int j = 0; j < dst.cols; j++)
		{
			dst(i, j)[0] = Y(i, j);
			dst(i, j)[1] = U(i, j);
			dst(i, j)[2] = V(i, j);
		}
	}

	return dst;
}

long int writeInFile(vector<vector<vector<pair<int, int>>>> data1)
{
	long int size = 0;
	ofstream myfile;
	myfile.open("compressed.txt");

	for (vector<vector<pair<int, int>>> data : data1)
	{
		int n = data.size();
		myfile << n << " "; //number of blocks
		size += 1;

		for (vector<pair<int, int>> block_data : data)
		{
			int m = block_data.size();
			myfile << m << " ";
			size += 1;

			for (pair<int, int> pair_data : block_data)
			{
				myfile << pair_data.first << " " << pair_data.second << " ";
				size += 1;
			}
		}
	}

	myfile.close();

	return size;
}

vector<vector<vector<pair<int, int>>>> readFromFile()
{
	ifstream myfile;
	myfile.open("compressed.txt");
	int n = 1, m = 1, x = 0, y = 0;
	vector<vector<vector<pair<int, int>>>> image_data;
	bool ok = false;
	while (!myfile.eof())
	{
		vector<vector<pair<int, int>>> data;

		myfile >> n; //number of blocks

		for (int i = 0; i < n; i++)
		{
			myfile >> m;
			vector<pair<int, int>> block_data;

			for (int j = 0; j < m; j++)
			{
				myfile >> x >> y;
				block_data.push_back(std::make_pair(x, y));
			}

			data.push_back(block_data);
		}

		image_data.push_back(data);
	}

	return image_data;
}


int main()
{
	Mat_<Vec3b> src = imread("Images/baboon_bmp.bmp", CV_LOAD_IMAGE_COLOR);

	uchar luminanceValues[65] = { 16, 11, 10, 16, 24, 40, 51, 61,
								12, 12, 14, 19, 26, 58, 60, 55,
								14, 13, 16, 24, 40, 57, 69, 56,
								14, 17, 22, 29, 51, 87, 80, 62,
								18, 22, 37, 56, 68, 109, 103, 77,
								24, 35, 55, 64, 81, 104, 113, 92,
								49, 64, 78, 87, 103, 121, 120, 101,
								72, 92, 95, 98, 112, 100, 103, 99 };

	uchar chrominanceValues[65] = { 17, 18, 24, 47, 99, 99, 99, 99,
									18, 21, 26, 66, 99, 99, 99, 99,
									24, 26, 56, 99, 99, 99, 99, 99,
									47, 66, 99, 99, 99, 99, 99, 99,
									99, 99, 99, 99, 99, 99, 99, 99,
									99, 99, 99, 99, 99, 99, 99, 99,
									99, 99, 99, 99, 99, 99, 99, 99,
									99, 99, 99, 99, 99, 99, 99, 99 };

	Mat_<uchar> luminanceQuantization(8, 8, luminanceValues);
	Mat_<uchar> chrominanceQuantization(8, 8, chrominanceValues);

	Mat img_out;

	cvtColor(src, img_out, CV_BGR2YCrCb);

	imshow("initial", src);
	imshow("yuv", img_out);

	vector<Mat_<uchar>> spl = split(img_out);

	imshow("initial Y", spl[0]);
	imshow("initial U", spl[1]);
	imshow("initial V", spl[2]);

	vector<vector<pair<int, int>>> compressOfU = compress(spl[1], chrominanceQuantization);

	vector<vector<pair<int, int>>> compressOfV = compress(spl[2], chrominanceQuantization);

	vector<vector<pair<int, int>>> compressOfY = compress(spl[0], luminanceQuantization);

	vector<vector<vector<pair<int, int>>>> data1;

	data1.push_back(compressOfY);
	data1.push_back(compressOfU);
	data1.push_back(compressOfV);

	long int fileSize = writeInFile(data1);
	long int fileSizeInKb = (fileSize * 4) / 1024;

	vector<vector<vector<pair<int, int>>>> image_data = readFromFile();

	Mat_<uchar> decompressOfU = decompress(spl[1].rows, spl[1].cols, image_data[1], chrominanceQuantization);
	Mat_<uchar> decompressOfV = decompress(spl[2].rows, spl[2].cols, image_data[2], chrominanceQuantization);
	Mat_<uchar> decompressOfY = decompress(spl[0].rows, spl[0].cols, image_data[0], luminanceQuantization);

	imshow("decompress of Y", decompressOfY);
	imshow("decompress of U", decompressOfU);
	imshow("decompress of V", decompressOfV);

	Mat_<Vec3b> decompressed = reassemble(decompressOfY, decompressOfU, decompressOfV);

	Mat_<Vec3b> result;
	//convert to RGB
	cvtColor(decompressed, result, CV_YCrCb2BGR);

	imshow("decompressed", result);

	imwrite("compressed.jpg", result);

	size_t originalSizeInKb = (src.total() * src.elemSize()) / 1024;

	cout << "original size " << originalSizeInKb << "\n";
	cout << "compressed size " << fileSizeInKb << "\n";

	float compression_ratio = (float)fileSizeInKb / originalSizeInKb;
	cout << "compression ratio " << compression_ratio;

	waitKey();

	return 0;
}