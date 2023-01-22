// 3次元WE-FDTD法+PML CUDA version
// 任意境界条件version
// 2017.4.06
// ver.0.12(PML0.1)
// Takao Tsuchiya, Doshisha Univ.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <mpi.h>
#include <sys/stat.h>

//#define M_PI 3.14159265358979
#define HBD_TO_HBU 1
#define HBU_TO_HBD 10

// モデルデータ
int iVox;					// モデルタイプ
char InputName[200] = {}, ModelName[200] = {};	// 入力ファイル名，モデル名
char VoxName[200] = {};			// ボクセルファイル名
char ObsName[10][200] = {};			// 観測点ファイル名
int Scheme;					// 手法
int3 Ndiv, Ndim;			// x, y, z方向分割数
int Nzorg;					// 元のz方向分割数
int Nreg, Nt;				// 領域分割数, 計算ステップ数
int Nobs, Nfo, Nsrc, Nfobs[10];			// 観測点数
int Mobs = 150000;
struct Pnt{					// 観測点情報
	int x, y, z;
	float t, p;
};
Pnt Src;					// 音源座標
Pnt* obs; 					// 観測点座標
int3* hobs, ocnt[10];		 	// 観測点数, 観測点中心座標
float* drv;				 	// 音源波形
float cfl, fs, dl, dt;			// CFL, サンプリング周波数, dl
float Ref[7], Ref0;			// 境界反射率
int wl;						// 入力波長
int iplane, ipn, iptime, iwave;		// 出力平面，平面位置，時間間隔，波形出力
int Nwave;					// 波形データ一括転送時間ステップ数
char WaveName[200] = {};			// 波形データファイル名
float* pp;					// 音圧分布図

// for GPU
int Ngpu;					// 使用GPU数
int GpuId = -1;				// 1GPUのときのGPU ID
int3 Block;					// Blockサイズ
int Bblock = 256;
int Boff;					// Block offset
float mem;

// for MPI
int Nnode, inode;			// ノード数，ランク
MPI_Status mpi_stat;		// MPIステータス

// for PML
int Nl;						// PML層数(0でMur))
int3 Nd;
float alpha, c1, c2, c3;
float c = 340;

#include "WE-FDTD.h"


int main(int argc, char* argv[]) 
{

	// モデルデータ読み込み
	ReadModel(argc, argv);

	// CE-FDTD設定値
	float4 d;
	
	if(cfl == 0) cfl = 1. / sqrt(3);
	
	d.x = cfl * cfl;
	d.w = 2. * (1. - 3. * cfl * cfl);
	dt = cfl * dl / c;
	fs = 1.0 / dt;

	mem = (double)(Ndiv.x)*Ndiv.y*Ndiv.z*8 + (double)Nwave*Nobs*28 + Nobs*12*7;
	
	if(Nl == 0)
		printf("SLF method + Mur\n");
	else{
		printf("SLF method + PML\n");
		printf(" Number of PML layer = %d\n", Nl);
		Nd.x = Nl * Ndim.y * Ndim.z;
		Nd.y = Ndiv.x * Nl * Ndim.z;
		Nd.z = Ndiv.x * Ndiv.y * Nl;
		c1 = c * c * dt / dl;
		c2 = dt / dl;
		c3 = dt / dl;
//			alpha = 0;
//		alpha = 100000 * dt;
//		alpha = 200000 * dt;
		mem += (double)(Nd.x+Nd.y+Nd.z)*48;
	}
	if(iVox > 0){
		printf(" Model file name: %s\n", VoxName);
		mem += (double)(Ndiv.x)*Ndiv.y*Ndiv.z;
	}
	else
		printf(" Rectangular model\n");
		
	printf(" CFL = %f, Fs = %f(Hz)\n", cfl, fs);
	printf(" dl = %f(m), dt = %e(s)\n", dl, dt);
	printf(" Nx = %d, Ny = %d, Nz = %d, Nt = %d\n", Ndiv.x, Ndiv.y, Ndiv.z, Nt);
	printf(" Size: %f x %f x %f(m) = %f(m^3), %f(s)\n", Ndiv.x*dl, Ndiv.y*dl, Ndiv.z*dl, 
		Ndiv.x*dl*Ndiv.y*dl*Ndiv.z*dl, Nt/fs);
	Ref0 = (cfl - 1.0) / (cfl + 1.0);


	// 1ブロックBlock_x*Block_yのスレッド数で並列計算
	int bdx = Ndiv.x / Block.x;
	int bdy = Ndiv.y / Block.y;
	int bdz = Ndiv.z / Block.z;
	printf(" Block: %d, %d, %d,", Block.x, Block.y, Block.z);
	printf(" Block Size: %d * %d * %d, Thread Size: %d\n", bdx, bdy, bdz, Block.x*Block.y*Block.z);

	// 観測波形用
	Nwave = 100;
	float* wave  = (float*) malloc(sizeof(float)*Nwave*Nobs);	// 観測点音圧波形
	float* hwave = (float*) malloc(sizeof(float)*Nwave*Nobs*7);	// 観測点音圧波形
	float* uu  = (float*) malloc(sizeof(float)*Nobs);			// 観測点音圧波形
	for(int i = 0; i < Nobs; i++)
		uu[i] = 0;
	for(int i = 0; i < Nwave*Nobs; i++)
		wave[i] = 0;
	for(int i = 0; i < Nwave*Nobs*7; i++)
		hwave[i] = 0;

	printf(" Obserbation points:\n");
	if(Nfo > 0){
		for(int ifo = 0; ifo < Nfo; ifo++){
			printf("  %s: no. points = %d\n", ObsName[ifo], Nfobs[ifo]);
		}
	}
	printf("  Total obs. points = %d\n", Nobs);
	
	// "data" ディレクトリ検査，作成
	int errd = 0;
	struct stat st;
	char DirName[200] = "data/";
	
	errd = stat(DirName, &st);
	if(errd == -1)
		errd = mkdir(DirName, 0777);

	strcpy(WaveName, DirName);
	if(argc == 1){
		strcat(WaveName, "wave");
//		strcpy(WaveName, "wave");
	}
	else{
		strcat(WaveName, "wave_");
//		strcpy(WaveName, "wave_");
		int len = strlen(InputName);
		strncat(WaveName, InputName, len-4);
	}
	
	FILE *fpo[10], *fp2;
	if(Nfo > 0){
		for(int ifo = 0; ifo < Nfo; ifo++){
			char TmpName[200] = {};
			strcpy(TmpName, WaveName);
			strcat(TmpName, "_");
			int len = strlen(ObsName[ifo]);
			strncat(TmpName, ObsName[ifo], len-4);
			sprintf(TmpName, "%s%d", TmpName, ifo);
			if(iwave == 0)
				strcat(TmpName, ".csv");
			else
				strcat(TmpName, ".bin");
			printf(" Ooutput (wave) file: %s\n", TmpName);
			if(iwave == 0){
				fpo[ifo] = fopen(TmpName,"w");
			}
			else{
				fpo[ifo] = fopen(TmpName,"wb");
				fwrite(&Nfobs[ifo], sizeof(int), 1, fpo[ifo]);
				fwrite(&Nt, sizeof(int), 1, fpo[ifo]);
			}
		}
	}
	else{
		if(iwave == 0)
			strcat(WaveName, ".csv");
		else
			strcat(WaveName, ".bin");
		printf(" Ooutput (wave) file: %s\n", WaveName);
		if(iwave == 0){
			fp2 = fopen(WaveName,"w");
		}
		else{
			fp2 = fopen(WaveName,"wb");
			fwrite(&Nobs, sizeof(int), 1, fp2);
			fwrite(&Nt, sizeof(int), 1, fp2);
		}
	}
	printf("\n");


	// 音圧分布図用配列
	if(iplane == 1)		// xy
		pp = (float*) malloc(sizeof(float)*Ndiv.x*Ndiv.y);
	if(iplane == 2)		// yz
		pp = (float*) malloc(sizeof(float)*Ndiv.y*Ndiv.z);
	if(iplane == 3)		// xz
		pp = (float*) malloc(sizeof(float)*Ndiv.x*Ndiv.z);


	int Num_gpu = 0, gpu_id;
	cudaGetDeviceCount(&Num_gpu);
	cudaEvent_t start,stop;
	cudaSetDevice(GpuId);	// "% num_gpus" allows more CPU threads than GPU devices
	cudaGetDevice(&gpu_id);
	printf(" %d GPUs found, No. %d device is used\n", Num_gpu, gpu_id);

	unsigned char* Vox;
	unsigned long long *Bid;
	unsigned short *Bnode;
	int Nbnd = 0;
	unsigned long long id, Nem;
	int nbx;
	float* pobs  = (float*) malloc(sizeof(float));			 	// 観測点音圧

	Nem = (unsigned long long)Ndiv.x * Ndiv.y * Ndiv.z;

	if(iVox > 0){
		Vox  = (unsigned char*)malloc(sizeof(unsigned char)*Nem);
		for(int k = 0; k < Ndiv.z; k++){
			for(int j = 0; j < Ndiv.y; j++){
				for(int i = 0; i < Ndiv.x; i++){
					id = (unsigned long long)Ndiv.x * Ndiv.y * k + Ndiv.x * j + i;
					Vox[id] = 0;
				}
			}
		}
		Nbnd = VoXRead(Vox, ModelName);
		nbx = Nbnd / Bblock + 1;
		Nbnd = Bblock * nbx;

		Bid = (unsigned long long*)malloc(sizeof(unsigned long long)*Nbnd);
		Bnode = (unsigned short*)malloc(sizeof(unsigned short)*Nbnd);
		for(int i = 0; i < Nbnd; i++){
			Bid[i] = 0;
			Bnode[i] = 0;
		}
		SetBoundary(Vox, Bid, Bnode, Nem);
		mem += Nbnd*10;
	}
	
//	float* mwave = (float*) malloc(sizeof(float)*Nwave*Nobs*7);	// 観測点音圧波形
//	for(int i = 0; i < Nwave*Nobs*7; i++)
//		mwave[i] = 0;

	cudaDeviceProp dev;
	cudaGetDeviceProperties(&dev, gpu_id);
	printf(" Global Memory Usage: %f (GB), Total Global Memory %f (GB)\n\n", mem/1024./1024./1024., 
		dev.totalGlobalMem/1024./1024./1024.);
	if(mem > dev.totalGlobalMem){
		printf(" Momory over!!\n");
		exit(1);
	}

	// デバイス上にメモリを確保する
	float *dp, *dpp, *tmp, *u, *pm, *dRef;
	float *dwave;
	unsigned char *dVox;
	unsigned long long *dBid;
	unsigned short *dBnode;
	int3 *dobs;				// 観測点座標
	
	cudaMalloc((void**) &dp,   sizeof(float)*Nem);			// 音圧
	cudaMalloc((void**) &dpp,  sizeof(float)*Nem);			// 1ステップ前音圧
	cudaMalloc((void**) &dRef, sizeof(float)*7);			// 反射係数
	if(iVox > 0){
		cudaMalloc((void**) &dVox, sizeof(unsigned char)*Nem);			// 形状データ
		cudaMalloc((void**) &dBid, sizeof(unsigned long long)*Nbnd);	// 境界条件idデータ
		cudaMalloc((void**) &dBnode, sizeof(unsigned short)*Nbnd);		// 境界反射係数データ
	}
	cudaMalloc((void**) &dobs, sizeof(int3)*Nobs*7);			// 観測波形用
	cudaMalloc((void**) &dwave, sizeof(float)*Nwave*Nobs*7);	// 観測波形用
	
	
	// デバイスメモリの初期化
	cudaMemset(dp,   0, sizeof(float)*Nem);
	cudaMemset(dpp,  0, sizeof(float)*Nem);
	cudaMemcpy(dRef, Ref, sizeof(float)*7, cudaMemcpyHostToDevice);
	cudaMemcpy(dobs, hobs, sizeof(int3)*Nobs*7, cudaMemcpyHostToDevice);
	cudaMemset(dwave,0, sizeof(float)*Nwave*Nobs*7);
	if(iVox > 0){
		cudaMemcpy(dVox, Vox, sizeof(unsigned char)*Nem, cudaMemcpyHostToDevice);
		cudaMemcpy(dBid, Bid, sizeof(unsigned long long)*Nbnd, cudaMemcpyHostToDevice);
		cudaMemcpy(dBnode, Bnode, sizeof(unsigned short)*Nbnd, cudaMemcpyHostToDevice);
	}
	if(Nl > 0){
		cudaMalloc((void**) &u,  sizeof(float)*6*(Nd.x+Nd.y+Nd.z));
		cudaMalloc((void**) &pm, sizeof(float)*6*(Nd.x+Nd.y+Nd.z));
		cudaMemset(u,  0, sizeof(float)*6*(Nd.x+Nd.y+Nd.z));
		cudaMemset(pm, 0, sizeof(float)*6*(Nd.x+Nd.y+Nd.z));
	}
	dim3 grid(bdx, bdy, bdz);
	dim3 threads(Block.x, Block.y, Block.z);
	dim3 gridb(nbx, 1, 1);
	dim3 threadsb(Bblock, 1, 1);
	dim3 gridw(1, 1, 1);
	dim3 threadsw(1, 1, 1);
	
	cudaError_t err = cudaGetLastError();
	printf("cuda debug:: line:%d rank:%d gpu:%d msg:%s\n", __LINE__, inode, gpu_id, cudaGetErrorString(err));

	float4 Driv;
	Driv.x = Driv.y = Driv.z = Driv.w = 0.0;
	float theta, phi;
	theta = M_PI / 180. * Src.t;	// 仰角
	phi   = M_PI / 180. * Src.p;	// 水平角
	int itt;
	
	// タイマーを作成して計測開始
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);

	// 時間ループ
	for(int it = 0; it < Nt; it++){

		itt = it % 100;
		if(itt == 0 && inode == 0){
			if(Ngpu == 1 || gpu_id == 0)
				printf("step: %d\n", it);
			err = cudaGetLastError();
			if(err != 0)
				printf("cuda debug:: line:%d rank:%d gpu:%d msg:%s\n", __LINE__, inode, gpu_id, cudaGetErrorString(err));
		}

		// 音源計算
		Driv.w = drv[it] / 2.0;
		if(theta == 0.0 && phi == 0.0){
			Driv.w = drv[it];
		}
		else{
			Driv.x += drv[it] * cfl / 2. * sin(theta) * cos(phi);
			Driv.y += drv[it] * cfl / 2. * sin(theta) * sin(phi);
			Driv.z += drv[it] * cfl / 2. * cos(theta);
		}

		// 領域計算
		if(iVox > 0){
			WE_Vox<<<grid, threads>>>(dp, dpp, Ndiv, Ndim, Src, Driv, cfl, d, dVox, Nl);
		}
		else{
			WE<<<grid, threads>>>(dp, dpp, Ndiv, Ndim, Src, Driv, cfl, d, Nl);
		}
//		cudaThreadSynchronize();
		tmp = dpp; dpp = dp; dp = tmp;
		
		WaveObss<<<gridw, threadsw>>>(dp, Ndiv, dobs, dwave, Nobs, Nwave, itt, Nl);
//			cudaThreadSynchronize();

		// 境界計算
		if(iVox == 0){
			if(Nl == 0){
				WE_Mur<<<grid, threads>>>(dp, dpp, Ndiv, Ndim, Ref0, Nl);
			}
			else{
				WE_press0<<<grid, threads>>>(dp, dpp, Ndiv, Ndim, cfl, Nl);
				PML_press<<<grid, threads>>>(dp, u, pm, Ndiv, Ndim, Nd, c1, alpha, Nl);
				PML_particle<<<grid, threads>>>(dp, u, Ndiv, Ndim, Nd, c2, alpha, Nl);
			}
		}
		else{
			CE_boundary_Plane<<<gridb, threadsb>>>(dp, dpp, dRef, dVox, dBid, dBnode, Ndiv, cfl, Nbnd);
//			cudaThreadSynchronize();
			CE_boundary_Edge<<<gridb, threadsb>>>(dp, dpp, dRef, dVox, dBid, dBnode, Ndiv, cfl, Nbnd);
//			cudaThreadSynchronize();
			CE_boundary_Corner<<<gridb, threadsb>>>(dp, dpp, dRef, dVox, dBid, dBnode, Ndiv, cfl, Nbnd);
//			cudaThreadSynchronize();
			if(Nl > 0){
				WE_press0<<<grid, threads>>>(dp, dpp, Ndiv, Ndim, cfl, Nl);
				PML_press<<<grid, threads>>>(dp, u, pm, Ndiv, Ndim, Nd, c1, alpha, Nl);
				PML_particle<<<grid, threads>>>(dp, u, Ndiv, Ndim, Nd, c2, alpha, Nl);
			}
//				WE_Mur<<<grid, threads>>>(dp, dpp, Ndiv, Ndim, Ref0, Nl);
		}
//		cudaThreadSynchronize();
		
		// 観測点音圧取得
		double ux, uy, uz, theta2, phi2;
		int noo = 0;
		if(itt == Nwave - 1){
			cudaMemcpy(hwave, dwave, sizeof(float)*Nwave*Nobs*7, cudaMemcpyDeviceToHost);
			for(int ii = 0; ii < Nwave; ii++){
				for(int io = 0; io < Nobs*7; io++){
					if(abs(hwave[ii*Nobs*7+io]) > 10.0 || isnan(hwave[ii*Nobs*7+io]) != 0){
						printf(" Diverged! %d: %e\n", it, hwave[ii*Nobs*7+io]);
						exit(1);
					}
				}
			}
					
			noo = 0;
			for(int ii = 0; ii < Nwave; ii++){
				for(int io = 0; io < Nobs; io++){
					theta2 = obs[io].t / 180. * M_PI;
					phi2   = obs[io].p / 180. * M_PI;
				
					if(theta2 == 0.0 && phi2 == 0.0){
						wave[ii*Nobs+io] = hwave[noo];
					}
					else{
						ux = (hwave[noo+1] - hwave[noo+2]) * cfl / 2.;
						uy = (hwave[noo+3] - hwave[noo+4]) * cfl / 2.;
						uz = (hwave[noo+5] - hwave[noo+6]) * cfl / 2.;
						uu[io] = uu[io] + ux * sin(theta2) * cos(phi2) 
									   + uy * sin(theta2) * sin(phi2) + uz * cos(theta2);
						wave[ii*Nobs+io] = (uu[io] + hwave[noo]) / 2.;
					}
					noo += 7;
				}
			}

			int iobs = 0;
			if(iwave == 0){
				if(Nfo > 0){
					for(int ifo = 0; ifo < Nfo; ifo++){
						for(int ii = 0; ii < Nwave; ii++){
							for(int io = 0; io < Nfobs[ifo]; io++){
								fprintf(fpo[ifo], "%e,", wave[ii*Nobs+io+iobs]);
							}
							fprintf(fpo[ifo], "\n");
						}
						iobs += Nfobs[ifo];
					}
				}
				else{
					for(int ii = 0; ii < Nwave; ii++){
						for(int io = 0; io < Nobs; io++){
							fprintf(fp2, "%e,", wave[ii*Nobs+io]);
						}
						fprintf(fp2, "\n");
					}
				}
			}
			else{
				if(Nfo > 0){
					iobs = 0;
					for(int ifo = 0; ifo < Nfo; ifo++){
						for(int ii = 0; ii < Nwave; ii++){
							fwrite(wave+ii*Nobs+iobs, sizeof(float), Nfobs[ifo], fpo[ifo]);
						}
						iobs += Nfobs[ifo];
					}
				}
				else
					fwrite(wave, sizeof(float), Nwave*Nobs, fp2);
			}
		}
//		cudaThreadSynchronize();

		// 音圧分布図保存
		if(ipn > 0 && it % iptime == 0){
			save_cross_section(dp, pp, it, DirName);
		}

	}
	
	//タイマーを停止しかかった時間を表示
	if(inode == 0 && (Nreg == 1 || gpu_id == 0)){
		float elapsed_time = 0.0;
		cudaEventRecord(stop,0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&elapsed_time, start, stop);
		printf("time: %f s\n", elapsed_time / 1000.);
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
	}
	
//	fclose(fp2);
    return 0;

}

