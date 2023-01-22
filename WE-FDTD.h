void ReadModel(int argc, char* argv[])
{
	char buf[200], SrcName[20];
	int dum;

	if(argc == 1)
		strcpy(InputName, "input.dat");
	else
		strcpy(InputName, argv[1]);

	FILE *fi  = fopen(InputName,"r");
	if(fi == NULL){
		printf("error:: No input file!\n");
		exit(1);
	}

	if(fgets(buf, sizeof(buf), fi) != NULL)
	sscanf(buf, "%d", &dum);							// 1行目はコメント
	
	if(fgets(buf, sizeof(buf), fi) != NULL)
	sscanf(buf, "%d %d", &iVox, &Scheme);		// モデルタイプ(0: 矩形, 1: 任意), 手法(1:SLF,6:IWB), PML層数(0:Mur,>0:PML)
	if(Scheme != 1){
		printf("error:: Invalid scheme (SLF only)!\n");
		exit(1);
	}

	if(iVox > 0){										// 任意ボクセルモデル
		if(fgets(buf, sizeof(buf), fi) != NULL)
		sscanf(buf, "%s", ModelName);					// ボクセルデータファイル名
		int dum;
		
		sprintf(VoxName, "%s.vox", ModelName);	// ボクセルファイル
		FILE *fim  = fopen(VoxName, "rb");
		if(fim == NULL){
			if(inode == 0) printf("error:: No VOX file!\n");
			exit(1);
		}
		if(fgets(buf, sizeof(buf), fim) != NULL)
		sscanf(buf, "%d %d %d", &Ndiv.x, &Ndiv.y, &Ndiv.z);	// 分割数だけ先に読み込む
		if(fgets(buf, sizeof(buf), fim) != NULL)
		sscanf(buf, "%d %d %d", &dum, &dum, &dum);		// ダミー
		if(fgets(buf, sizeof(buf), fim) != NULL)
		sscanf(buf, "%d %d %d", &Block.x, &Block.y, &Block.z);
		fclose(fim);
	}
	else{												// 矩形モデル
		if(fgets(buf, sizeof(buf), fi) != NULL)
		sscanf(buf, "%d %d %d", &Ndiv.x, &Ndiv.y, &Ndiv.z);		// 分割数
		Block.x = 4;
		Block.y = 4;
		Block.z = 4;
	}
	Boff = Block.z / 2;

	for(int ir = 0; ir < 7; ir++)
		Ref[ir] = 0;
	if(fgets(buf, sizeof(buf), fi) != NULL)						// 境界0-境界5までの反射率，境界6は無反射
	sscanf(buf, "%f %f %f %f %f %f", &Ref[0], &Ref[1], &Ref[2], &Ref[3], &Ref[4], &Ref[5]);

	if(fgets(buf, sizeof(buf), fi) != NULL)
	sscanf(buf, "%d", &Nt);								// 計算ステップ数
	
	if(fgets(buf, sizeof(buf), fi) != NULL)
	sscanf(buf, "%d %d %d %f %f", &Src.x, &Src.y, &Src.z, &Src.t, &Src.p);	// 音源位置，方向読込
	if(Src.x > Ndiv.x || Src.y > Ndiv.y || Src.z > Ndiv.z){
		printf("error:: Invalid source location (%d, %d, %d)\n", Src.x, Src.y, Src.z);
		exit(1);
	}
	
	if(fgets(buf, sizeof(buf), fi) != NULL)
	sscanf(buf, "%d %s", &wl, SrcName);				// 音源波長，音源ファイル(0: source.csv, 1: インパルス, -1: 微分インパルス)

	drv = (float*) malloc(sizeof(float)*Nt); 			// 音源波形
	for(int it = 0; it < Nt; it++){
		drv[it] = 0.;
	}
	if(wl == 0){										// 音源ファイルから波形入力
		FILE *fi2  = fopen(SrcName, "r");
		float s = 0;
		for(int it = 0; it < Nt; it++){
			if(fgets(buf, sizeof(buf), fi2) != NULL){
				sscanf(buf, "%f", &drv[it]);
				s += drv[it];
			}
		}
		fclose(fi2);
		float ss = 0;
		for(int it = 0; it < Nt; it++){
			drv[it] -= s / Nt;
			ss += drv[it];
			drv[it] = ss;
		}
	}
	else if(wl == 1){		// インパルス
		drv[0] = 1.;
	}
	else if(wl == -1){		// 微分インパルス
		drv[0] = 1.;
		drv[1] = -1.;
	}
	else{					// 正弦波1波パルス
		int nw = wl;
		if(wl > Nt) nw = Nt-1;
		for(int it = 0; it <= nw; it++){
			drv[it] = sin(2 * 3.1415926 * it / wl);
//			printf("%f\n", drv[it]);
		}
	}

	if(fgets(buf, sizeof(buf), fi) != NULL)
	sscanf(buf, "%f %f", &cfl, &dl);					// CFL数，dl
	
	int ifo = 0, iob, Nob, ico;
	int io = 0, ioo = 0;
	char TempName[200] = {};
	obs  = (Pnt*) malloc(sizeof(Pnt)*Mobs); 			// 観測点
	hobs = (int3*) malloc(sizeof(int3)*Mobs*7); 		// 観測点座標
	float ox, oy, oz;

	Nobs = 0;
	if(fgets(buf, sizeof(buf), fi) != NULL)
	sscanf(buf, "%d %d %d", &iob, &Nob, &ico);			// 観測ファイル数，観測点数
	if(iob > 0){
		Nfo = Nob;
		for(ifo = 0; ifo < Nfo; ifo++){
			ioo = 0;
			if(ico == 0){
				if(fgets(buf, sizeof(buf), fi) != NULL)
				sscanf(buf, "%s %f %f %f", TempName, &ox, &oy, &oz);		// 
				ocnt[ifo].x = (int)(ox / dl);
				ocnt[ifo].y = (int)(oy / dl);
				ocnt[ifo].z = (int)(oz / dl);
//				printf("%d %d %d\n", ocnt[ifo].x,ocnt[ifo].y,ocnt[ifo].z);
			}
			else{
				if(fgets(buf, sizeof(buf), fi) != NULL)
				sscanf(buf, "%s %d %d %d", TempName, &ocnt[ifo].x, &ocnt[ifo].y, &ocnt[ifo].z);		// 
			}
			for(int ic = 0; ic < 200; ic++)
				ObsName[ifo][ic] = TempName[ic];
			sprintf(ObsName[ifo], "%s.csv", ObsName[ifo]);	// 観測ファイル名
//			printf(" %s: %d %d %d\n", ObsName[ifo], ocnt[ifo].x, ocnt[ifo].y, ocnt[ifo].z);

			FILE *fobs  = fopen(ObsName[ifo], "r");
			if(fobs == NULL){
				printf("error:: No observe file! %s\n", ObsName[ifo]);
				exit(1);
			}
			while (!feof(fobs) && io < Mobs){
				if(ico == 0){
					// 観測点のx,y,z座標，仰角θ, 水平角φ
					if(fgets(buf, sizeof(buf), fobs) != NULL)
					sscanf(buf, "%f,%f,%f,%f,%f", &ox, &oy, &oz, &obs[io].t, &obs[io].p);
					obs[io].x = (int)(ox / dl) + ocnt[ifo].x;
					obs[io].y = (int)(oy / dl) + ocnt[ifo].y;
					obs[io].z = (int)(oz / dl) + ocnt[ifo].z;
				}
				else{
					// 観測点のx,y,z座標，仰角θ, 水平角φ
					if(fgets(buf, sizeof(buf), fobs) != NULL)
					sscanf(buf, "%d,%d,%d,%f,%f", &obs[io].x, &obs[io].y, &obs[io].z, &obs[io].t, &obs[io].p);
					obs[io].x += ocnt[ifo].x;
					obs[io].y += ocnt[ifo].y;
					obs[io].z += ocnt[ifo].z;
				}

				if(obs[io].x < 1 || obs[io].x > Ndiv.x || obs[io].y < 1 || obs[io].y > Ndiv.y 
					|| obs[io].z < 1 || obs[io].z > Ndiv.z){
					printf("error:: Illegal observation point in %s, No.%d: %d %d %d\n", ObsName[ifo], ioo+1, obs[io].x, obs[io].y, obs[io].z);
					MPI_Abort( MPI_COMM_WORLD, 1 );
					exit(1);
				}
				for(int ii = 0; ii < 7; ii++){
					hobs[io*7+ii].x = obs[io].x;
					hobs[io*7+ii].y = obs[io].y;
					hobs[io*7+ii].z = obs[io].z;
				}
				hobs[io*7+1].x = obs[io].x+1;
				hobs[io*7+2].x = obs[io].x-1;
				hobs[io*7+3].y = obs[io].y+1;
				hobs[io*7+4].y = obs[io].y-1;
				hobs[io*7+5].z = obs[io].z+1;
				hobs[io*7+6].z = obs[io].z-1;
				io++;
				ioo++;
			}
			if(io > Mobs){
				printf("error:: No of observation points exceed %d!\n", Mobs);
				MPI_Abort( MPI_COMM_WORLD, 1 );
				exit(1);
			}
			--io;
			Nfobs[ifo] = ioo - 1;
			fclose(fobs);
			Nobs += Nfobs[ifo];
		}
	}
//		for(int i = 0; i < Nfo; i++)
//			printf("%d\n", Nfobs[i]);
//		printf("%d\n", Nobs);
	else{
		Nobs = Nob;
		for(int io = 0; io < Nobs; io++){
			if(ico == 0){
				if(fgets(buf, sizeof(buf), fi) != NULL)
				sscanf(buf, "%f %f %f %f %f", &ox, &oy, &oz, &obs[io].t, &obs[io].p);
				obs[io].x = (int)(ox / dl);
				obs[io].y = (int)(oy / dl);
				obs[io].z = (int)(oz / dl);
			}
			else{
				if(fgets(buf, sizeof(buf), fi) != NULL)
				// 観測点のx,y,z座標，仰角θ, 水平角φ
				sscanf(buf, "%d %d %d %f %f", &obs[io].x, &obs[io].y, &obs[io].z, &obs[io].t, &obs[io].p);
			}
			if(obs[io].x < 1 || obs[io].x > Ndiv.x || obs[io].y < 1 || obs[io].y > Ndiv.y 
				|| obs[io].z < 1 || obs[io].z > Ndiv.z){
				printf("error:: Illegal observation point, No.%d: %d %d %d\n", io+1, obs[io].x, obs[io].y, obs[io].z);
				MPI_Abort( MPI_COMM_WORLD, 1 );
				exit(1);
			}
			for(int ii = 0; ii < 7; ii++){
				hobs[io*7+ii].x = obs[io].x;
				hobs[io*7+ii].y = obs[io].y;
				hobs[io*7+ii].z = obs[io].z;
			}
			hobs[io*7+1].x = obs[io].x+1;
			hobs[io*7+2].x = obs[io].x-1;
			hobs[io*7+3].y = obs[io].y+1;
			hobs[io*7+4].y = obs[io].y-1;
			hobs[io*7+5].z = obs[io].z+1;
			hobs[io*7+6].z = obs[io].z-1;
		}
	}

	if(fgets(buf, sizeof(buf), fi) != NULL)
	sscanf(buf, "%d %d %d %d", &iplane, &ipn, &iptime, &iwave);		// コンター面，波形出力
	
	if(fgets(buf, sizeof(buf), fi) != NULL)
	sscanf(buf, "%d %d", &Ngpu, &GpuId);				// GPU数/ノード, 1GPUのときの使用GPU ID
	if(fgets(buf, sizeof(buf), fi) != NULL)
	sscanf(buf, "%d %f", &Nl, &alpha);					// PML層数，PMLアルファ
	
	Ndim.x = Ndiv.x - 2 * Nl;
	Ndim.y = Ndiv.y - 2 * Nl;
	Ndim.z = Ndiv.z - 2 * Nl;

	Ngpu = 1;
	Nnode = 1;
	Nreg = Nnode * Ngpu;
		
	// 分割チェック
	int Fx, Fy, Fz;
	Fx = Ndiv.x % Block.x;
	Fy = Ndiv.y % Block.y;
	Fz = Ndiv.z % Block.z;
	if(Fx + Fy + Fz > 0){
		if( inode == 0 ){
			printf("error:: Invalid number of divisions.\n" );
			if(Fx > 0) printf("Nx: %d is invalid\n", Ndiv.x);
			if(Fy > 0) printf("Ny: %d is invalid\n", Ndiv.y);
			if(Fz > 0) printf("Nz: %d is invalid\n", Ndiv.z);
		}
		exit(1);
	}

}


int VoXRead(unsigned char* Vox, char* ModelName)
{
	unsigned long long id;
	int Nbnd = 0, dum;
	char bufm[200];
		
	FILE *fim  = fopen(VoxName, "rb");		// ボクセルデータ
	if(fim == NULL){
		printf("error:: No VOX file!\n");
		exit(1);
	}

	// ボクセルデータ読み込み
	printf(" Reading voxel data: %s\n", VoxName);
	short* nn = (short*)malloc(sizeof(short)*Ndiv.x);
	unsigned char* in = (unsigned char*)malloc(sizeof(unsigned char)*Ndiv.x);
	short lnn;
	size_t siz;

	for(int ii = 0; ii < 5; ii++){
		if(fgets(bufm, sizeof(bufm), fim) != NULL);
		sscanf(bufm, "%d", &dum);
	}

	int ix;
	for(int k = 0; k < Ndiv.z; k++){
		for(int j = 0; j < Ndiv.y; j++){
			siz = fread(&lnn, sizeof(short), 1, fim);
			siz = fread(nn, sizeof(short), lnn, fim);
			siz = fread(in, sizeof(unsigned char), lnn, fim);
			if(siz == 0){
				printf("error:: Illegal boundary data!!\n");
				exit(1);
			}
			ix = 0;
			for(int i = 0; i < lnn; i++){
				for(int ii = 0; ii < nn[i]; ii++){
					id = (unsigned long long)Ndiv.x * Ndiv.y * k + Ndiv.x * j + ix;
					Vox[id] = in[i];
					ix++;
					if(Vox[id] > 1) Nbnd++;		// 境界条件データ数
				}
			}
		}
	}
	fclose(fim);

	if(Nbnd < 1){
		printf("error:: No boundary data!!\n");
		exit(1);
	}
	return Nbnd;
}


void SetBoundary(unsigned char *Vox, unsigned long long *Bid, unsigned short *Bnode, unsigned long long Nem)
{
	unsigned long long id = 0;
	int node, n[30];
	int Nx, Ny, Nz;
	int bx, by, bz, bb, bid = 0;

	Nx = Ndiv.x; Ny = Ndiv.y; Nz = Ndiv.z;
	for(int k = 1; k < Nz-1; k++){
		for(int j = 1; j < Ny-1; j++){
			for(int i = 1; i < Nx-1; i++){
				id = (unsigned long long)Nx * Ny * k + Nx * j + i;
				node = 0;
				if(Vox[id] > 1){
					n[0] = Vox[id-1         ];
					n[1] = Vox[id  +Nx      ];
					n[2] = Vox[id+1         ];
					n[3] = Vox[id  -Nx      ];
					n[4] = Vox[id     -Nx*Ny];
					n[5] = Vox[id     +Nx*Ny];
					n[6] = Vox[id-1-Nx      ];
					n[7] = Vox[id-1+Nx      ];
					n[8] = Vox[id+1+Nx      ];
					n[9] = Vox[id+1-Nx      ];
					n[10] = Vox[id  -Nx-Nx*Ny];
					n[11] = Vox[id  -Nx+Nx*Ny];
					n[12] = Vox[id  +Nx+Nx*Ny];
					n[13] = Vox[id  +Nx-Nx*Ny];
					n[14] = Vox[id-1   -Nx*Ny];
					n[15] = Vox[id-1   +Nx*Ny];
					n[16] = Vox[id+1   +Nx*Ny];
					n[17] = Vox[id+1   -Nx*Ny];
					n[18] = Vox[id-1-Nx-Nx*Ny];
					n[19] = Vox[id-1+Nx-Nx*Ny];
					n[20] = Vox[id+1+Nx-Nx*Ny];
					n[21] = Vox[id+1-Nx-Nx*Ny];
					n[22] = Vox[id-1-Nx+Nx*Ny];
					n[23] = Vox[id-1+Nx+Nx*Ny];
					n[24] = Vox[id+1+Nx+Nx*Ny];
					n[25] = Vox[id+1-Nx+Nx*Ny];

					bx = by = bz = bb = 0;

					if(n[0] == 1) {bx = 1; bb = 1;}
					if(n[1] == 1) {by = 2; bb = 1;}
					if(n[2] == 1) {bx = 2; bb = 1;}
					if(n[3] == 1) {by = 1; bb = 1;}
					if(n[4] == 1) {bz = 1; bb = 1;}
					if(n[5] == 1) {bz = 2; bb = 1;}

					if(bb == 0){
						bb = 0;
							 if(n[6]  == 1){bx = 1; by = 1;}
						else if(n[7]  == 1){bx = 1; by = 2;}
						else if(n[8]  == 1){bx = 2; by = 2;}
						else if(n[9]  == 1){bx = 2; by = 1;}
						else if(n[10] == 1){by = 1; bz = 1;}
						else if(n[11] == 1){by = 1; bz = 2;}
						else if(n[12] == 1){by = 2; bz = 2;}
						else if(n[13] == 1){by = 2; bz = 1;}
						else if(n[14] == 1){bx = 1; bz = 1;}
						else if(n[15] == 1){bx = 1; bz = 2;}
						else if(n[16] == 1){bx = 2; bz = 2;}
						else if(n[17] == 1){bx = 2; bz = 1;}

						else if(n[18] == 1){bx = 1; by = 1; bz = 1;}
						else if(n[19] == 1){bx = 1; by = 2; bz = 1;}
						else if(n[20] == 1){bx = 2; by = 2; bz = 1;}
						else if(n[21] == 1){bx = 2; by = 1; bz = 1;}
						else if(n[22] == 1){bx = 1; by = 1; bz = 2;}
						else if(n[23] == 1){bx = 1; by = 2; bz = 2;}
						else if(n[24] == 1){bx = 2; by = 2; bz = 2;}
						else if(n[25] == 1){bx = 2; by = 1; bz = 2;}
					}

					if((bx+by+bz) == 0){
						node = 0;
					}
					else{
						node = 27 * bb + 9 * bz + 3 * by + bx;
					}
					if(node > 0){
						id = (unsigned long long)k * Nx * Ny + Nx * j + i;
						Bid[bid] = id;
						Bnode[bid] = (unsigned short)node;
						++bid;
						if(id > Nem) printf("%d %d %d: %lld\n", i, j, k, id);
					}
				}
			}
		}
	}
	return;
}


__global__ void WE_Vox(float* dp, float* dpp, int3 Ndiv, int3 Ndim, Pnt Src, float4 Driv, float cfl, float4 d, 
		unsigned char* dVox, int Nl){
	
	const unsigned int tx  = threadIdx.x;
	const unsigned int ty  = threadIdx.y;
	const unsigned int tz  = threadIdx.z;
	const unsigned int bdx = blockDim.x;
	const unsigned int bdy = blockDim.y;
	const unsigned int bdz = blockDim.z;
	const unsigned int bx  = blockIdx.x;
	const unsigned int by  = blockIdx.y;
	const unsigned int bz  = blockIdx.z;

	const unsigned int ix = bx * bdx + tx;
	const unsigned int iy = by * bdy + ty;
	const unsigned int iz = bz * bdz + tz;
	unsigned long long id = Ndiv.x * Ndiv.y * iz + Ndiv.x * iy + ix;
//	unsigned long long iv = Ndiv.x * Ndiv.y * (iz-Nl) + Ndiv.x * (iy-Nl) + (ix-Nl);
	float4 dr;

	if(ix > Nl && ix < Ndim.x+Nl-1 && iy > Nl && iy < Ndim.y+Nl-1 && iz > Nl && iz < Ndim.z+Nl-1){
		if(dVox[id] == 1){
			dpp[id] = (dp[id-1] + dp[id+1] + dp[id-Ndiv.x] + dp[id+Ndiv.x] + dp[id-Ndiv.x*Ndiv.y] + dp[id+Ndiv.x*Ndiv.y]) * d.x
					+ d.w * dp[id] - dpp[id];

			dr.w = Driv.w - Driv.x - Driv.y - Driv.z;
			dr.x = Driv.x;
			dr.y = Driv.y;
			dr.z = Driv.z;

		if(ix == Src.x && iy == Src.y && iz == Src.z)
			dpp[id] += dr.w;
		if(ix == Src.x+1 && iy == Src.y && iz == Src.z)
			dpp[id] += dr.w;
		if(ix == Src.x && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.w;
		if(ix == Src.x && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.w;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.w;
		if(ix == Src.x+1 && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.w;
		if(ix == Src.x && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.w;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.w;

		if(ix == Src.x+1 && iy == Src.y && iz == Src.z)
			dpp[id] += dr.x;
		if(ix == Src.x+2 && iy == Src.y && iz == Src.z)
			dpp[id] += dr.x;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.x;
		if(ix == Src.x+1 && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.x;
		if(ix == Src.x+2 && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.x;
		if(ix == Src.x+2 && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.x;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.x;
		if(ix == Src.x+2 && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.x;

		if(ix == Src.x && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.y;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.y;
		if(ix == Src.x && iy == Src.y+2 && iz == Src.z)
			dpp[id] += dr.y;
		if(ix == Src.x && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.y;
		if(ix == Src.x+1 && iy == Src.y+2 && iz == Src.z)
			dpp[id] += dr.y;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.y;
		if(ix == Src.x && iy == Src.y+2 && iz == Src.z+1)
			dpp[id] += dr.y;
		if(ix == Src.x+1 && iy == Src.y+2 && iz == Src.z+1)
			dpp[id] += dr.y;
			
		if(ix == Src.x && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.z;
		if(ix == Src.x+1 && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.z;
		if(ix == Src.x && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.z;
		if(ix == Src.x && iy == Src.y && iz == Src.z+2)
			dpp[id] += dr.z;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.z;
		if(ix == Src.x+1 && iy == Src.y && iz == Src.z+2)
			dpp[id] += dr.z;
		if(ix == Src.x && iy == Src.y+1 && iz == Src.z+2)
			dpp[id] += dr.z;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z+2)
			dpp[id] += dr.z;
		}
	}
 	__syncthreads();
 	
}


__global__ void WE(float* dp, float* dpp, int3 Ndiv, int3 Ndim, Pnt Src, float4 Driv, float cfl, float4 d, int Nl){
	
	const unsigned int tx  = threadIdx.x;
	const unsigned int ty  = threadIdx.y;
	const unsigned int tz  = threadIdx.z;
	const unsigned int bdx = blockDim.x;
	const unsigned int bdy = blockDim.y;
	const unsigned int bdz = blockDim.z;
	const unsigned int bx  = blockIdx.x;
	const unsigned int by  = blockIdx.y;
	const unsigned int bz  = blockIdx.z;

	const unsigned int ix = bx * bdx + tx;
	const unsigned int iy = by * bdy + ty;
	const unsigned int iz = bz * bdz + tz;
	unsigned long long id = Ndiv.x * Ndiv.y * iz + Ndiv.x * iy + ix;

	float4 dr;
	
	if(ix > Nl && ix < Ndim.x+Nl-1 && iy > Nl && iy < Ndim.y+Nl-1 && iz > Nl && iz < Ndim.z+Nl-1){
		dpp[id] = (dp[id-1] + dp[id+1] + dp[id-Ndiv.x] + dp[id+Ndiv.x] + dp[id-Ndiv.x*Ndiv.y] + dp[id+Ndiv.x*Ndiv.y]) * d.x
				+ d.w * dp[id] - dpp[id];

		dr.w = Driv.w - Driv.x - Driv.y - Driv.z;
		dr.x = Driv.x;
		dr.y = Driv.y;
		dr.z = Driv.z;

		if(ix == Src.x && iy == Src.y && iz == Src.z)
			dpp[id] += dr.w;
		if(ix == Src.x+1 && iy == Src.y && iz == Src.z)
			dpp[id] += dr.w;
		if(ix == Src.x && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.w;
		if(ix == Src.x && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.w;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.w;
		if(ix == Src.x+1 && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.w;
		if(ix == Src.x && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.w;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.w;

		if(ix == Src.x+1 && iy == Src.y && iz == Src.z)
			dpp[id] += dr.x;
		if(ix == Src.x+2 && iy == Src.y && iz == Src.z)
			dpp[id] += dr.x;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.x;
		if(ix == Src.x+1 && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.x;
		if(ix == Src.x+2 && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.x;
		if(ix == Src.x+2 && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.x;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.x;
		if(ix == Src.x+2 && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.x;

		if(ix == Src.x && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.y;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z)
			dpp[id] += dr.y;
		if(ix == Src.x && iy == Src.y+2 && iz == Src.z)
			dpp[id] += dr.y;
		if(ix == Src.x && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.y;
		if(ix == Src.x+1 && iy == Src.y+2 && iz == Src.z)
			dpp[id] += dr.y;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.y;
		if(ix == Src.x && iy == Src.y+2 && iz == Src.z+1)
			dpp[id] += dr.y;
		if(ix == Src.x+1 && iy == Src.y+2 && iz == Src.z+1)
			dpp[id] += dr.y;
			
		if(ix == Src.x && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.z;
		if(ix == Src.x+1 && iy == Src.y && iz == Src.z+1)
			dpp[id] += dr.z;
		if(ix == Src.x && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.z;
		if(ix == Src.x && iy == Src.y && iz == Src.z+2)
			dpp[id] += dr.z;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z+1)
			dpp[id] += dr.z;
		if(ix == Src.x+1 && iy == Src.y && iz == Src.z+2)
			dpp[id] += dr.z;
		if(ix == Src.x && iy == Src.y+1 && iz == Src.z+2)
			dpp[id] += dr.z;
		if(ix == Src.x+1 && iy == Src.y+1 && iz == Src.z+2)
			dpp[id] += dr.z;
	}
 	__syncthreads();
 	
}


__global__ void WaveObss(float* dp, int3 Ndiv, int3* dobs, float* dwave, int Nobs, int Nwave, int itt, int Nl){
	
	unsigned long long id;

	for(int io = 0; io < Nobs*7; io++){
		id = (dobs[io].z % Ndiv.z) * Ndiv.x * Ndiv.y + dobs[io].y * Ndiv.x + dobs[io].x;
			dwave[itt*Nobs*7+io] = dp[id];
	}
 	
}


__global__ void CE_boundary_Plane(float* dp, float* dpp, float* dRef, unsigned char* dVox, 
	unsigned long long *dBid, unsigned short *dBnode, int3 Ndiv, float cfl, int Nbnd){
	
	const unsigned int tx  = threadIdx.x;
	const unsigned int bdx = blockDim.x;
	const unsigned int bx  = blockIdx.x;
	const unsigned int tid = bx * bdx + tx;
	const unsigned int Nx = Ndiv.x;
	const unsigned int Ny = Ndiv.y;

	unsigned long long id;
	int n, br, brz, bry, brx, nref, ir;
	double ref;
	unsigned char v, bit;

	if(tid < Nbnd){
		id = dBid[tid];
		if(id > 0){
			n = dBnode[tid];
			v = dVox[id];
			if(n > 0 && v > 1){
				ref = 0;
				bit = 2;
				nref = 0;
				for(ir = 0; ir < 7; ir++){
					ref += (v & bit)/bit * dRef[ir];
					nref += (v & bit)/bit;
					bit = bit * 2;
				}
				if(nref > 0) ref = ref / nref;

				ref = ((1.0 + ref) * cfl - (1.0 - ref)) / ((1.0 + ref) * cfl + (1.0 - ref));
				br  = (int)(n / 27);
				brz = (int)((n - br*27) / 9);
				bry = (int)((n - br*27 - brz*9) / 3);
				brx = n % 3;

				if(brx == 1 && bry == 0 && brz == 0){
					dp[id] = dpp[id-1] + ref * (dp[id-1] - dpp[id]);
				}
				if(brx == 2 && bry == 0 && brz == 0){
					dp[id] = dpp[id+1] + ref * (dp[id+1] - dpp[id]);
				}
				if(brx == 0 && bry == 1 && brz == 0){
					dp[id] = dpp[id-Nx] + ref * (dp[id-Nx] - dpp[id]);
				}
				if(brx == 0 && bry == 2 && brz == 0){
					dp[id] = dpp[id+Nx] + ref * (dp[id+Nx] - dpp[id]);
				}
				if(brx == 0 && bry == 0 && brz == 1){
					dp[id] = dpp[id-Nx*Ny] + ref * (dp[id-Nx*Ny] - dpp[id]);
				}
				if(brx == 0 && bry == 0 && brz == 2){
					dp[id] = dpp[id+Nx*Ny] + ref * (dp[id+Nx*Ny] - dpp[id]);
				}
			}
		}
	}
	__syncthreads();
}


__global__ void CE_boundary_Edge(float* dp, float* dpp, float* dRef, unsigned char* dVox, 
	unsigned long long *dBid, unsigned short *dBnode, int3 Ndiv, float cfl, int Nbnd){
	
	const unsigned int tx  = threadIdx.x;
	const unsigned int bdx = blockDim.x;
	const unsigned int bx  = blockIdx.x;
	const unsigned int tid = bx * bdx + tx;
	const unsigned int Nx = Ndiv.x;
	const unsigned int Ny = Ndiv.y;

	unsigned long long id;
	int n, br, brz, bry, brx, jx, jy, jz, nref, ir;
	double ref, ref1, ref2, ref3, px, py, pz, pt;
	double pbx, pby, pbz;
	unsigned char v, bit;

	if(tid < Nbnd){
		id = dBid[tid];
		if(id > 0){
			n = dBnode[tid];
			v = dVox[id];
			if(n > 0 && v > 1){
				ref = 0;
				bit = 2;
				nref = 0;
				for(ir = 0; ir < 7; ir++){
					ref += (v & bit)/bit * dRef[ir];
					nref += (v & bit)/bit;
					bit = bit * 2;
				}
				if(nref > 0) ref = ref / nref;

				ref1 = (1.0 - ref) / (sqrt(2.0) * (1.0 + ref) * cfl + (1.0 - ref));
				ref2 = (1.0 + ref) * cfl /sqrt(2.0) / (sqrt(2.0) * (1.0 + ref) * cfl + (1.0 - ref));
				ref3 = ((1.0 + ref) * cfl - (1.0 - ref)) / ((1.0 + ref) * cfl + (1.0 - ref));
				br  = (int)(n / 27);
				brz = (int)((n - br*27) / 9);
				bry = (int)((n - br*27 - brz*9) / 3);
				brx = n % 3;

				if(br == 0){
					if(brx == 1 && bry == 1 && brz == 0){
						jx = -1;
						jy = -Nx;
						px = dp[id+jy]+dpp[id+jy]+dpp[id]-(dp[id+jx]+dp[id+jx+jy]+dpp[id+jx]+dpp[id+jx+jy]);
						py = dp[id+jx]+dpp[id+jx]+dpp[id]-(dp[id+jy]+dp[id+jx+jy]+dpp[id+jy]+dpp[id+jx+jy]);
						pt = dpp[id]+dpp[id+jx]+dpp[id+jy]+dpp[id+jx+jy]-(dp[id+jx]+dp[id+jy]+dp[id+jx+jy]);
						dp[id] = ref1 * pt - ref2 *(px + py);
					}
					if(brx == 2 && bry == 1 && brz == 0){
						jx =  1;
						jy = -Nx;
						px = dp[id+jy]+dpp[id+jy]+dpp[id]-(dp[id+jx]+dp[id+jx+jy]+dpp[id+jx]+dpp[id+jx+jy]);
						py = dp[id+jx]+dpp[id+jx]+dpp[id]-(dp[id+jy]+dp[id+jx+jy]+dpp[id+jy]+dpp[id+jx+jy]);
						pt = dpp[id]+dpp[id+jx]+dpp[id+jy]+dpp[id+jx+jy]-(dp[id+jx]+dp[id+jy]+dp[id+jx+jy]);
						dp[id] = ref1 * pt - ref2 *(px + py);
					}
					if(brx == 1 && bry == 2 && brz == 0){
						jx = -1;
						jy =  Nx;
						px = dp[id+jy]+dpp[id+jy]+dpp[id]-(dp[id+jx]+dp[id+jx+jy]+dpp[id+jx]+dpp[id+jx+jy]);
						py = dp[id+jx]+dpp[id+jx]+dpp[id]-(dp[id+jy]+dp[id+jx+jy]+dpp[id+jy]+dpp[id+jx+jy]);
						pt = dpp[id]+dpp[id+jx]+dpp[id+jy]+dpp[id+jx+jy]-(dp[id+jx]+dp[id+jy]+dp[id+jx+jy]);
						dp[id] = ref1 * pt - ref2 *(px + py);
					}
					if(brx == 2 && bry == 2 && brz == 0){
						jx = 1;
						jy = Nx;
						px = dp[id+jy]+dpp[id+jy]+dpp[id]-(dp[id+jx]+dp[id+jx+jy]+dpp[id+jx]+dpp[id+jx+jy]);
						py = dp[id+jx]+dpp[id+jx]+dpp[id]-(dp[id+jy]+dp[id+jx+jy]+dpp[id+jy]+dpp[id+jx+jy]);
						pt = dpp[id]+dpp[id+jx]+dpp[id+jy]+dpp[id+jx+jy]-(dp[id+jx]+dp[id+jy]+dp[id+jx+jy]);
						dp[id] = ref1 * pt - ref2 *(px + py);
					}

					if(brx == 0 && bry == 1 && brz == 1){
						jy = -Nx;
						jz = -Nx*Ny;
						py = dp[id+jz]+dpp[id+jz]+dpp[id]-(dp[id+jy]+dp[id+jy+jz]+dpp[id+jy]+dpp[id+jy+jz]);
						pz = dp[id+jy]+dpp[id+jy]+dpp[id]-(dp[id+jz]+dp[id+jy+jz]+dpp[id+jz]+dpp[id+jy+jz]);
						pt = dpp[id]+dpp[id+jy]+dpp[id+jz]+dpp[id+jy+jz]-(dp[id+jy]+dp[id+jz]+dp[id+jy+jz]);
						dp[id] = ref1 * pt - ref2 *(py + pz);
					}
					if(brx == 0 && bry == 2 && brz == 1){
						jy =  Nx;
						jz = -Nx*Ny;
						py = dp[id+jz]+dpp[id+jz]+dpp[id]-(dp[id+jy]+dp[id+jy+jz]+dpp[id+jy]+dpp[id+jy+jz]);
						pz = dp[id+jy]+dpp[id+jy]+dpp[id]-(dp[id+jz]+dp[id+jy+jz]+dpp[id+jz]+dpp[id+jy+jz]);
						pt = dpp[id]+dpp[id+jy]+dpp[id+jz]+dpp[id+jy+jz]-(dp[id+jy]+dp[id+jz]+dp[id+jy+jz]);
						dp[id] = ref1 * pt - ref2 *(py + pz);
					}
					if(brx == 0 && bry == 1 && brz == 2){
						jy = -Nx;
						jz =  Nx*Ny;
						py = dp[id+jz]+dpp[id+jz]+dpp[id]-(dp[id+jy]+dp[id+jy+jz]+dpp[id+jy]+dpp[id+jy+jz]);
						pz = dp[id+jy]+dpp[id+jy]+dpp[id]-(dp[id+jz]+dp[id+jy+jz]+dpp[id+jz]+dpp[id+jy+jz]);
						pt = dpp[id]+dpp[id+jy]+dpp[id+jz]+dpp[id+jy+jz]-(dp[id+jy]+dp[id+jz]+dp[id+jy+jz]);
						dp[id] = ref1 * pt - ref2 *(py + pz);
					}
					if(brx == 0 && bry == 2 && brz == 2){
						jy = Nx;
						jz = Nx*Ny;
						py = dp[id+jz]+dpp[id+jz]+dpp[id]-(dp[id+jy]+dp[id+jy+jz]+dpp[id+jy]+dpp[id+jy+jz]);
						pz = dp[id+jy]+dpp[id+jy]+dpp[id]-(dp[id+jz]+dp[id+jy+jz]+dpp[id+jz]+dpp[id+jy+jz]);
						pt = dpp[id]+dpp[id+jy]+dpp[id+jz]+dpp[id+jy+jz]-(dp[id+jy]+dp[id+jz]+dp[id+jy+jz]);
						dp[id] = ref1 * pt - ref2 *(py + pz);
					}

					if(brx == 1 && bry == 0 && brz == 1){
						jx = -1;
						jz = -Nx*Ny;
						px = dp[id+jz]+dpp[id+jz]+dpp[id]-(dp[id+jx]+dp[id+jx+jz]+dpp[id+jx]+dpp[id+jx+jz]);
						pz = dp[id+jx]+dpp[id+jx]+dpp[id]-(dp[id+jz]+dp[id+jx+jz]+dpp[id+jz]+dpp[id+jx+jz]);
						pt = dpp[id]+dpp[id+jx]+dpp[id+jz]+dpp[id+jx+jz]-(dp[id+jx]+dp[id+jz]+dp[id+jx+jz]);
						dp[id] = ref1 * pt - ref2 *(px + pz);
					}
					if(brx == 2 && bry == 0 && brz == 1){
						jx =  1;
						jz = -Nx*Ny;
						px = dp[id+jz]+dpp[id+jz]+dpp[id]-(dp[id+jx]+dp[id+jx+jz]+dpp[id+jx]+dpp[id+jx+jz]);
						pz = dp[id+jx]+dpp[id+jx]+dpp[id]-(dp[id+jz]+dp[id+jx+jz]+dpp[id+jz]+dpp[id+jx+jz]);
						pt = dpp[id]+dpp[id+jx]+dpp[id+jz]+dpp[id+jx+jz]-(dp[id+jx]+dp[id+jz]+dp[id+jx+jz]);
						dp[id] = ref1 * pt - ref2 *(px + pz);
					}
					if(brx == 1 && bry == 0 && brz == 2){
						jx = -1;
						jz =  Nx*Ny;
						px = dp[id+jz]+dpp[id+jz]+dpp[id]-(dp[id+jx]+dp[id+jx+jz]+dpp[id+jx]+dpp[id+jx+jz]);
						pz = dp[id+jx]+dpp[id+jx]+dpp[id]-(dp[id+jz]+dp[id+jx+jz]+dpp[id+jz]+dpp[id+jx+jz]);
						pt = dpp[id]+dpp[id+jx]+dpp[id+jz]+dpp[id+jx+jz]-(dp[id+jx]+dp[id+jz]+dp[id+jx+jz]);
						dp[id] = ref1 * pt - ref2 *(px + pz);
					}
					if(brx == 2 && bry == 0 && brz == 2){
						jx = 1;
						jz = Nx*Ny;
						px = dp[id+jz]+dpp[id+jz]+dpp[id]-(dp[id+jx]+dp[id+jx+jz]+dpp[id+jx]+dpp[id+jx+jz]);
						pz = dp[id+jx]+dpp[id+jx]+dpp[id]-(dp[id+jz]+dp[id+jx+jz]+dpp[id+jz]+dpp[id+jx+jz]);
						pt = dpp[id]+dpp[id+jx]+dpp[id+jz]+dpp[id+jx+jz]-(dp[id+jx]+dp[id+jz]+dp[id+jx+jz]);
						dp[id] = ref1 * pt - ref2 *(px + pz);
					}
				}
				else{
					if(brx == 1 && bry == 1 && brz == 0){
						pbx = dpp[id-1] + ref3 * (dp[id-1] - dpp[id]);
						pby = dpp[id-Nx] + ref3 * (dp[id-Nx] - dpp[id]);
						dp[id] = (pbx + pby) / 2.0;
					}
					if(brx == 2 && bry == 1 && brz == 0){
						pbx = dpp[id+1] + ref3 * (dp[id+1] - dpp[id]);
						pby = dpp[id-Nx] + ref3 * (dp[id-Nx] - dpp[id]);
						dp[id] = (pbx + pby) / 2.0;
					}
					if(brx == 1 && bry == 2 && brz == 0){
						pbx = dpp[id-1] + ref3 * (dp[id-1] - dpp[id]);
						pby = dpp[id+Nx] + ref3 * (dp[id+Nx] - dpp[id]);
						dp[id] = (pbx + pby) / 2.0;
					}
					if(brx == 2 && bry == 2 && brz == 0){
						pbx = dpp[id+1] + ref3 * (dp[id+1] - dpp[id]);
						pby = dpp[id+Nx] + ref3 * (dp[id+Nx] - dpp[id]);
						dp[id] = (pbx + pby) / 2.0;
					}

					if(brx == 0 && bry == 1 && brz == 1){
						pby = dpp[id-Nx] + ref3 * (dp[id-Nx] - dpp[id]);
						pbz = dpp[id-Nx*Ny] + ref3 * (dp[id-Nx*Ny] - dpp[id]);
						dp[id] = (pby + pbz) / 2.0;
					}
					if(brx == 0 && bry == 2 && brz == 1){
						pby = dpp[id+Nx] + ref3 * (dp[id+Nx] - dpp[id]);
						pbz = dpp[id-Nx*Ny] + ref3 * (dp[id-Nx*Ny] - dpp[id]);
						dp[id] = (pby + pbz) / 2.0;
					}
					if(brx == 0 && bry == 1 && brz == 2){
						pby = dpp[id-Nx] + ref3 * (dp[id-Nx] - dpp[id]);
						pbz = dpp[id+Nx*Ny] + ref3 * (dp[id+Nx*Ny] - dpp[id]);
						dp[id] = (pby + pbz) / 2.0;
					}
					if(brx == 0 && bry == 2 && brz == 2){
						pby = dpp[id+Nx] + ref3 * (dp[id+Nx] - dpp[id]);
						pbz = dpp[id+Nx*Ny] + ref3 * (dp[id+Nx*Ny] - dpp[id]);
						dp[id] = (pby + pbz) / 2.0;
					}

					if(brx == 1 && bry == 0 && brz == 1){
						pbx = dpp[id-1] + ref3 * (dp[id-1] - dpp[id]);
						pbz = dpp[id-Nx*Ny] + ref3 * (dp[id-Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pbz) / 2.0;
					}
					if(brx == 2 && bry == 0 && brz == 1){
						pbx = dpp[id+1] + ref3 * (dp[id+1] - dpp[id]);
						pbz = dpp[id-Nx*Ny] + ref3 * (dp[id-Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pbz) / 2.0;
					}
					if(brx == 1 && bry == 0 && brz == 2){
						pbx = dpp[id-1] + ref3 * (dp[id-1] - dpp[id]);
						pbz = dpp[id+Nx*Ny] + ref3 * (dp[id+Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pbz) / 2.0;
					}
					if(brx == 2 && bry == 0 && brz == 2){
						pbx = dpp[id+1] + ref3 * (dp[id+1] - dpp[id]);
						pbz = dpp[id+Nx*Ny] + ref3 * (dp[id+Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pbz) / 2.0;
					}
				}
			}
		}
	}
 	__syncthreads();
}


__global__ void CE_boundary_Corner(float* dp, float* dpp, float* dRef, unsigned char* dVox, 
	unsigned long long *dBid, unsigned short *dBnode, int3 Ndiv, float cfl, int Nbnd){

	const unsigned int tx  = threadIdx.x;
	const unsigned int bdx = blockDim.x;
	const unsigned int bx  = blockIdx.x;
	const unsigned int tid = bx * bdx + tx;
	const unsigned int Nx = Ndiv.x;
	const unsigned int Ny = Ndiv.y;

	unsigned long long id;
	int n, br, brz, bry, brx, jx, jy, jz, nref, ir;
	double px, py, pz, pt, ref, ref1, ref2, ref3;
	double pbx = 0, pby = 0, pbz = 0;
	unsigned char v, bit;

	if(tid < Nbnd){
		id = dBid[tid];
		if(id > 0){
			n = dBnode[tid];
			v = dVox[id];
			if(n > 0 && v > 1){
				ref = 0;
				bit = 2;
				nref = 0;
				for(ir = 0; ir < 7; ir++){
					ref += (v & bit)/bit * dRef[ir];
					nref += (v & bit)/bit;
					bit = bit * 2;
				}
				if(nref > 0) ref = ref / nref;

				ref1 = (1.0 - ref) / (sqrt(3.0)*(1.0 + ref) * cfl + (1.0 - ref));
				ref2 = (1.0 + ref) * cfl / sqrt(3.0) / (sqrt(3.0)*(1.0 + ref) * cfl + (1.0 - ref));
				ref3 = ((1.0 + ref) * cfl - (1.0 - ref)) / ((1.0 + ref) * cfl + (1.0 - ref));
				br  = (int)(n / 27);
				brz = (int)((n - br*27) / 9);
				bry = (int)((n - br*27 - brz*9) / 3);
				brx = n % 3;

				if(br == 0){
					if(brx == 1 && bry == 1 && brz == 1){
						jx = -1;
						jy = -Nx;
						jz = -Nx * Ny;
						px = dp[id+jy+jz]+dp[id+jz]+dp[id+jy]+dpp[id+jy+jz]+dpp[id+jz]+dpp[id+jy]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jx]+dpp[id+jx+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jx]);
						py = dp[id+jx+jz]+dp[id+jz]+dp[id+jx]+dpp[id+jx+jz]+dpp[id+jz]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jy]+dp[id+jy]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jy]+dpp[id+jy]);
						pz = dp[id+jx+jy]+dp[id+jy]+dp[id+jx]+dpp[id+jx+jy]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jz]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jz]);
						pt = dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jz]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jz]+dp[id+jy]+dp[id+jx]);
						dp[id] = ref1 * pt - ref2 *(px + py + pz);
					}
					if(brx == 2 && bry == 1 && brz == 1){
						jx =  1;
						jy = -Nx;
						jz = -Nx * Ny;
						px = dp[id+jy+jz]+dp[id+jz]+dp[id+jy]+dpp[id+jy+jz]+dpp[id+jz]+dpp[id+jy]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jx]+dpp[id+jx+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jx]);
						py = dp[id+jx+jz]+dp[id+jz]+dp[id+jx]+dpp[id+jx+jz]+dpp[id+jz]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jy]+dp[id+jy]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jy]+dpp[id+jy]);
						pz = dp[id+jx+jy]+dp[id+jy]+dp[id+jx]+dpp[id+jx+jy]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jz]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jz]);
						pt = dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jz]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jz]+dp[id+jy]+dp[id+jx]);
						dp[id] = ref1 * pt - ref2 *(px + py + pz);
					}
					if(brx == 1 && bry == 2 && brz == 1){
						jx = -1;
						jy =  Nx;
						jz = -Nx * Ny;
						px = dp[id+jy+jz]+dp[id+jz]+dp[id+jy]+dpp[id+jy+jz]+dpp[id+jz]+dpp[id+jy]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jx]+dpp[id+jx+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jx]);
						py = dp[id+jx+jz]+dp[id+jz]+dp[id+jx]+dpp[id+jx+jz]+dpp[id+jz]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jy]+dp[id+jy]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jy]+dpp[id+jy]);
						pz = dp[id+jx+jy]+dp[id+jy]+dp[id+jx]+dpp[id+jx+jy]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jz]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jz]);
						pt = dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jz]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jz]+dp[id+jy]+dp[id+jx]);
						dp[id] = ref1 * pt - ref2 *(px + py + pz);
					}
					if(brx == 1 && bry == 1 && brz == 2){
						jx = -1;
						jy = -Nx;
						jz =  Nx * Ny;
						px = dp[id+jy+jz]+dp[id+jz]+dp[id+jy]+dpp[id+jy+jz]+dpp[id+jz]+dpp[id+jy]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jx]+dpp[id+jx+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jx]);
						py = dp[id+jx+jz]+dp[id+jz]+dp[id+jx]+dpp[id+jx+jz]+dpp[id+jz]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jy]+dp[id+jy]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jy]+dpp[id+jy]);
						pz = dp[id+jx+jy]+dp[id+jy]+dp[id+jx]+dpp[id+jx+jy]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jz]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jz]);
						pt = dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jz]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jz]+dp[id+jy]+dp[id+jx]);
						dp[id] = ref1 * pt - ref2 *(px + py + pz);
					}
					if(brx == 2 && bry == 2 && brz == 1){
						jx =  1;
						jy =  Nx;
						jz = -Nx * Ny;
						px = dp[id+jy+jz]+dp[id+jz]+dp[id+jy]+dpp[id+jy+jz]+dpp[id+jz]+dpp[id+jy]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jx]+dpp[id+jx+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jx]);
						py = dp[id+jx+jz]+dp[id+jz]+dp[id+jx]+dpp[id+jx+jz]+dpp[id+jz]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jy]+dp[id+jy]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jy]+dpp[id+jy]);
						pz = dp[id+jx+jy]+dp[id+jy]+dp[id+jx]+dpp[id+jx+jy]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jz]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jz]);
						pt = dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jz]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jz]+dp[id+jy]+dp[id+jx]);
						dp[id] = ref1 * pt - ref2 *(px + py + pz);
					}
					if(brx == 2 && bry == 1 && brz == 2){
						jx =  1;
						jy = -Nx;
						jz =  Nx * Ny;
						px = dp[id+jy+jz]+dp[id+jz]+dp[id+jy]+dpp[id+jy+jz]+dpp[id+jz]+dpp[id+jy]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jx]+dpp[id+jx+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jx]);
						py = dp[id+jx+jz]+dp[id+jz]+dp[id+jx]+dpp[id+jx+jz]+dpp[id+jz]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jy]+dp[id+jy]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jy]+dpp[id+jy]);
						pz = dp[id+jx+jy]+dp[id+jy]+dp[id+jx]+dpp[id+jx+jy]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jz]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jz]);
						pt = dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jz]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jz]+dp[id+jy]+dp[id+jx]);
						dp[id] = ref1 * pt - ref2 *(px + py + pz);
					}
					if(brx == 1 && bry == 2 && brz == 2){
						jx = -1;
						jy =  Nx;
						jz =  Nx * Ny;
						px = dp[id+jy+jz]+dp[id+jz]+dp[id+jy]+dpp[id+jy+jz]+dpp[id+jz]+dpp[id+jy]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jx]+dpp[id+jx+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jx]);
						py = dp[id+jx+jz]+dp[id+jz]+dp[id+jx]+dpp[id+jx+jz]+dpp[id+jz]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jy]+dp[id+jy]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jy]+dpp[id+jy]);
						pz = dp[id+jx+jy]+dp[id+jy]+dp[id+jx]+dpp[id+jx+jy]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jz]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jz]);
						pt = dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jz]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jz]+dp[id+jy]+dp[id+jx]);
						dp[id] = ref1 * pt - ref2 *(px + py + pz);
					}
					if(brx == 2 && bry == 2 && brz == 2){
						jx = 1;
						jy = Nx;
						jz = Nx * Ny;
						px = dp[id+jy+jz]+dp[id+jz]+dp[id+jy]+dpp[id+jy+jz]+dpp[id+jz]+dpp[id+jy]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jx]+dpp[id+jx+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jx]);
						py = dp[id+jx+jz]+dp[id+jz]+dp[id+jx]+dpp[id+jx+jz]+dpp[id+jz]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jy]+dp[id+jy]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jy]+dpp[id+jy]);
						pz = dp[id+jx+jy]+dp[id+jy]+dp[id+jx]+dpp[id+jx+jy]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jz]+dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jz]);
						pt = dpp[id+jx+jy+jz]+dpp[id+jy+jz]+dpp[id+jx+jz]+dpp[id+jx+jy]+dpp[id+jz]+dpp[id+jy]+dpp[id+jx]+dpp[id]
						   -(dp[id+jx+jy+jz]+dp[id+jy+jz]+dp[id+jx+jz]+dp[id+jx+jy]+dp[id+jz]+dp[id+jy]+dp[id+jx]);
						dp[id] = ref1 * pt - ref2 *(px + py + pz);
					}
				}
				else
				{
					if(brx == 1 && bry == 1 && brz == 1){
						pbx = dpp[id-1] + ref3 * (dp[id-1] - dpp[id]);
						pby = dpp[id-Nx] + ref3 * (dp[id-Nx] - dpp[id]);
						pbz = dpp[id-Nx*Ny] + ref3 * (dp[id-Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pby + pbz) / 3.0;
					}
					if(brx == 2 && bry == 1 && brz == 1){
						pbx = dpp[id+1] + ref3 * (dp[id+1] - dpp[id]);
						pby = dpp[id-Nx] + ref3 * (dp[id-Nx] - dpp[id]);
						pbz = dpp[id-Nx*Ny] + ref3 * (dp[id-Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pby + pbz) / 3.0;
					}
					if(brx == 1 && bry == 2 && brz == 1){
						pbx = dpp[id-1] + ref3 * (dp[id-1] - dpp[id]);
						pby = dpp[id+Nx] + ref3 * (dp[id+Nx] - dpp[id]);
						pbz = dpp[id-Nx*Ny] + ref3 * (dp[id-Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pby + pbz) / 3.0;
					}
					if(brx == 1 && bry == 1 && brz == 2){
						pbx = dpp[id-1] + ref3 * (dp[id-1] - dpp[id]);
						pby = dpp[id-Nx] + ref3 * (dp[id-Nx] - dpp[id]);
						pbz = dpp[id+Nx*Ny] + ref3 * (dp[id+Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pby + pbz) / 3.0;
					}
					if(brx == 2 && bry == 2 && brz == 1){
						pbx = dpp[id+1] + ref3 * (dp[id+1] - dpp[id]);
						pby = dpp[id+Nx] + ref3 * (dp[id+Nx] - dpp[id]);
						pbz = dpp[id-Nx*Ny] + ref3 * (dp[id-Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pby + pbz) / 3.0;
					}
					if(brx == 2 && bry == 1 && brz == 2){
						pbx = dpp[id+1] + ref3 * (dp[id+1] - dpp[id]);
						pby = dpp[id-Nx] + ref3 * (dp[id-Nx] - dpp[id]);
						pbz = dpp[id+Nx*Ny] + ref3 * (dp[id+Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pby + pbz) / 3.0;
					}
					if(brx == 1 && bry == 2 && brz == 2){
						pbx = dpp[id-1] + ref3 * (dp[id-1] - dpp[id]);
						pby = dpp[id+Nx] + ref3 * (dp[id+Nx] - dpp[id]);
						pbz = dpp[id+Nx*Ny] + ref3 * (dp[id+Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pby + pbz) / 3.0;
					}
					if(brx == 2 && bry == 2 && brz == 2){
						pbx = dpp[id+1] + ref3 * (dp[id+1] - dpp[id]);
						pby = dpp[id+Nx] + ref3 * (dp[id+Nx] - dpp[id]);
						pbz = dpp[id+Nx*Ny] + ref3 * (dp[id+Nx*Ny] - dpp[id]);
						dp[id] = (pbx + pby + pbz) / 3.0;
					}
				}
			}
		}
	}
 	__syncthreads();
}


__global__ void WE_Mur(float* dp, float* dpp, int3 Ndiv, int3 Ndim, float Ref, int Nl)
{
	const unsigned int tx  = threadIdx.x;
	const unsigned int ty  = threadIdx.y;
	const unsigned int tz  = threadIdx.z;
	const unsigned int bdx = blockDim.x;
	const unsigned int bdy = blockDim.y;
	const unsigned int bdz = blockDim.z;
	const unsigned int bx  = blockIdx.x;
	const unsigned int by  = blockIdx.y;
	const unsigned int bz  = blockIdx.z;

	const unsigned int ix = bx * bdx + tx;
	const unsigned int iy = by * bdy + ty;
	const unsigned int iz = bz * bdz + tz;
	unsigned long long id = Ndiv.x * Ndiv.y * iz + Ndiv.x * iy + ix;

	if(ix == Nl){
		dp[id] = dpp[id+1] + Ref * dp[id+1] - Ref * dpp[id];
	}
	if(ix == Ndim.x+Nl-1){
		dp[id] = dpp[id-1] + Ref * dp[id-1] - Ref * dpp[id];
	}
	if(iy == Nl){
		dp[id] = dpp[id+Ndim.x] + Ref * dp[id+Ndim.x] - Ref * dpp[id];
	}
	if(iy == Ndim.y+Nl-1){
		dp[id] = dpp[id-Ndim.x] + Ref * dp[id-Ndim.x] - Ref * dpp[id];
	}
	if(iz == Nl){
		dp[id] = dpp[id+Ndim.x*Ndim.y] + Ref * dp[id+Ndim.x*Ndim.y] - Ref * dpp[id];
	}
	if(iz == Ndim.z+Nl-1){
		dp[id] = dpp[id-Ndim.x*Ndim.y] + Ref * dp[id-Ndim.x*Ndim.y] - Ref * dpp[id];
	}
 	__syncthreads();
}


__global__ void WE_press0(float* dp, float* dpp, int3 Ndiv, int3 Ndim, float cfl, int Nl){

	const unsigned int tx  = threadIdx.x;
	const unsigned int ty  = threadIdx.y;
	const unsigned int tz  = threadIdx.z;
	const unsigned int bdx = blockDim.x;
	const unsigned int bdy = blockDim.y;
	const unsigned int bdz = blockDim.z;
	const unsigned int bx  = blockIdx.x;
	const unsigned int by  = blockIdx.y;
	const unsigned int bz  = blockIdx.z;

	const unsigned int ix = bx * bdx + tx;
	const unsigned int iy = by * bdy + ty;
	const unsigned int iz = bz * bdz + tz;
	const unsigned int Nx = Ndim.x;
	const unsigned int Ny = Ndim.y;
	const unsigned int Nz = Ndim.z;
	const unsigned int Mx = Ndiv.x;
	const unsigned int My = Ndiv.y;
	const unsigned long long id = iz * Mx*My + iy * Mx + ix;
	
	if(ix == Nl && iy > Nl && iy < Ny+Nl-1 && iz > Nl && iz < Nz+Nl-1){
		dp[id] = 2*dpp[id] - dp[id] + cfl*cfl * (dpp[id+1]+dpp[id-1] + dpp[id-Mx]+dpp[id+Mx] + dpp[id-Mx*My]+dpp[id+Mx*My] - 6*dpp[id]);
	}
	if(ix == Nx+Nl-1 && iy > Nl && iy < Ny+Nl-1 && iz > Nl && iz < Nz+Nl-1){
		dp[id] = 2*dpp[id] - dp[id] + cfl*cfl * (dpp[id+1]+dpp[id-1] + dpp[id-Mx]+dpp[id+Mx] + dpp[id-Mx*My]+dpp[id+Mx*My] - 6*dpp[id]);
	}
	if(ix >= Nl && ix <= Nx+Nl-1 && iy == Nl && iz > Nl && iz < Nz+Nl-1){
		dp[id] = 2*dpp[id] - dp[id] + cfl*cfl * (dpp[id+1]+dpp[id-1] + dpp[id-Mx]+dpp[id+Mx] + dpp[id-Mx*My]+dpp[id+Mx*My] - 6*dpp[id]);
	}
	if(ix >= Nl && ix <= Nx+Nl-1 && iy == Ny+Nl-1 && iz > Nl && iz < Nz+Nl-1){
		dp[id] = 2*dpp[id] - dp[id] + cfl*cfl * (dpp[id+1]+dpp[id-1] + dpp[id-Mx]+dpp[id+Mx] + dpp[id-Mx*My]+dpp[id+Mx*My] - 6*dpp[id]);
	}
	if(ix >= Nl && ix <= Nx+Nl-1 && iy >= Nl && iy <= Ny+Nl-1 && iz == Nl){
		dp[id] = 2*dpp[id] - dp[id] + cfl*cfl * (dpp[id+1]+dpp[id-1] + dpp[id-Mx]+dpp[id+Mx] + dpp[id-Mx*My]+dpp[id+Mx*My] - 6*dpp[id]);
	}
	if(ix >= Nl && ix <= Nx+Nl-1 && iy >= Nl && iy <= Ny+Nl-1 && iz == Nz+Nl-1){
		dp[id] = 2*dpp[id] - dp[id] + cfl*cfl * (dpp[id+1]+dpp[id-1] + dpp[id-Mx]+dpp[id+Mx] + dpp[id-Mx*My]+dpp[id+Mx*My] - 6*dpp[id]);
	}

}


__global__ void PML_press(float* dp, float* u, float* pm, int3 Ndiv, int3 Ndim, int3 Nd, float c1, float alpha, int Nl){

	const unsigned int tx  = threadIdx.x;
	const unsigned int ty  = threadIdx.y;
	const unsigned int tz  = threadIdx.z;
	const unsigned int bdx = blockDim.x;
	const unsigned int bdy = blockDim.y;
	const unsigned int bdz = blockDim.z;
	const unsigned int bx  = blockIdx.x;
	const unsigned int by  = blockIdx.y;
	const unsigned int bz  = blockIdx.z;

	const unsigned int ix = bx * bdx + tx;
	const unsigned int iy = by * bdy + ty;
	const unsigned int iz = bz * bdz + tz;
	const unsigned int Nx = Ndim.x;
	const unsigned int Ny = Ndim.y;
	const unsigned int Nz = Ndim.z;
	const unsigned int Mx = Ndiv.x;
	const unsigned int My = Ndiv.y;
	const unsigned int Mz = Ndiv.z;
	const unsigned int Mxy = Ndiv.x*Ndiv.y;
	const unsigned int Ndx = 3*Nd.x;
	const unsigned int Ndy = 3*Nd.y;
	const unsigned int Ndz = 3*Nd.z;
	const unsigned int ip = iz * Mxy + iy * Mx + ix;
	long long im;
	float rx, ry, rz;


// -x
	if(ix < Nl && iy >= Nl && iy < Ny+Nl && iz >= Nl && iz < Nz+Nl){
		im = 3 * ((iz-Nl)*Nl*Ny + (iy-Nl) * Nl + ix);
		rx = 0; ry = 0; rz = 0;
		rx = alpha * pow((float)(Nl-ix-0.5)/Nl, 2) / 2;
		if(ix > 0)
			pm[im]   = pm[im]   * (1 - rx) / (1 + rx) - c1 * (u[im]   - u[im-3]) / (1 + rx);
		if(iy == Nl){
			pm[im+1] = pm[im+1] * (1 - ry) / (1 + ry) - c1 * (u[im+1] - u[2*Ndx+3*(ix+(Nl-1)*Mx+(iz-Nl)*Nl*Mx)+1]) / (1 + ry);
		}
		else if(iy == Ny+Nl-1){
			pm[im+1] = pm[im+1] * (1 - ry) / (1 + ry) - c1 * (u[2*Ndx+Ndy+3*(ix+(iz-Nl)*Nl*Mx)+1] - u[im-3*Nl+1]) / (1 + ry);
		}
		else{
			pm[im+1] = pm[im+1] * (1 - ry) / (1 + ry) - c1 * (u[im+1] - u[im-3*Nl+1]) / (1 + ry);
		}
		if(iz == Nl)
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[im+2] - u[2*Ndx+2*Ndy+3*(ix+iy*Mx)+3*(Nl-1)*Mxy+2]) / (1 + rz);
		else if(iz == Nz+Nl-1)
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[2*Ndx+2*Ndy+Ndz+3*(ix+iy*Mx)+2] - u[im-3*Nl*Ny+2]) / (1 + rz);
		else
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[im+2] - u[im-3*Nl*Ny+2]) / (1 + rz);
		dp[ip] = pm[im] + pm[im+1] + pm[im+2];
	}

// +x
	if(ix >= Nx+Nl && iy >= Nl && iy < Ny+Nl && iz >= Nl && iz < Nz+Nl){
		im = 3 * ((iz-Nl)*Nl*Ny + (iy-Nl)*Nl + ix-(Nx+Nl)) + Ndx;
		rx = 0; ry = 0; rz = 0;
		rx = alpha * pow((float)(ix-(Nx+Nl)+0.5)/Nl, 2) / 2;
		if(ix < Mx-1)
			pm[im]   = pm[im]   * (1 - rx) / (1 + rx) - c1 * (u[im+3] - u[im]) / (1 + rx);
		if(iy == Nl){
			pm[im+1] = pm[im+1] * (1 - ry) / (1 + ry) - c1 * (u[im+1] - u[2*Ndx+3*(ix+(Nl-1)*Mx+(iz-Nl)*Nl*Mx)+1]) / (1 + ry);
		}
		else if(iy == Ny+Nl-1){
			pm[im+1] = pm[im+1] * (1 - ry) / (1 + ry) - c1 * (u[2*Ndx+Ndy+3*(ix+(iz-Nl)*Nl*Mx)+1] - u[im-3*Nl+1]) / (1 + ry);
		}
		else{
			pm[im+1] = pm[im+1] * (1 - ry) / (1 + ry) - c1 * (u[im+1] - u[im-3*Nl+1]) / (1 + ry);
		}
		if(iz == Nl)
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[im+2] - u[2*Ndx+2*Ndy+3*(ix+iy*Mx)+3*(Nl-1)*Mxy+2]) / (1 + rz);
		else if(iz == Nz+Nl-1)
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[2*Ndx+2*Ndy+Ndz+3*(ix+iy*Mx)+2] - u[im-3*Nl*Ny+2]) / (1 + rz);
		else
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[im+2] - u[im-3*Nl*Ny+2]) / (1 + rz);
		dp[ip] = pm[im] + pm[im+1] + pm[im+2];
	}

// -y
	if(iy < Nl && iz >= Nl && iz < Nz+Nl){
		im = 3 * ((iz-Nl)*Nl*Mx + iy*Mx + ix) + 2*Ndx;
		rx = 0; ry = 0; rz = 0;
		if(ix <     Nl) rx = alpha * pow((float)(Nl-ix-0.5)/Nl, 2) / 2;
		if(ix >= Nx+Nl) rx = alpha * pow((float)(ix-Nx-Nl+0.5)/Nl, 2) / 2;
		ry = alpha * pow((float)(Nl-iy-0.5)/Nl, 2) / 2;
		if(ix > 0 && ix < Mx-1)
			pm[im]   = pm[im]   * (1 - rx) / (1 + rx) - c1 * (u[im] - u[im-3]) / (1 + rx);
		if(iy > 0)
			pm[im+1] = pm[im+1] * (1 - ry) / (1 + ry) - c1 * (u[im+1] - u[im-3*Mx+1]) / (1 + ry);
		if(iz == Nl)
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[im+2] - u[2*Ndx+2*Ndy+3*(ix+iy*Mx)+3*(Nl-1)*Mxy+2]) / (1 + rz);
		else if(iz == Nz+Nl-1)
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[2*Ndx+2*Ndy+Ndz+3*(ix+iy*Mx)+2] - u[im-3*Nl*Mx+2]) / (1 + rz);
		else
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[im+2] - u[im-3*Nl*Mx+2]) / (1 + rz);
		dp[ip] = pm[im] + pm[im+1] + pm[im+2];
	}

// +y
	if(iy >= Ny+Nl && iz >= Nl && iz < Nz+Nl){
		im = 3 * ((iz-Nl)*Nl*Mx + (iy-(Ny+Nl))*Mx + ix) + 2*Ndx + Ndy;
		rx = 0; ry = 0; rz = 0;
		if(ix <     Nl) rx = alpha * pow((float)(Nl-ix-0.5)/Nl, 2) / 2;
		if(ix >= Nx+Nl) rx = alpha * pow((float)(ix-Nx-Nl+0.5)/Nl, 2) / 2;
		ry = alpha * pow((float)(iy-(Ny+Nl)+0.5)/Nl, 2) / 2;
		if(ix > 0 && ix < Mx-1)
			pm[im]   = pm[im]   * (1 - rx) / (1 + rx) - c1 * (u[im] - u[im-3]) / (1 + rx);
		if(iy < My-1)
			pm[im+1] = pm[im+1] * (1 - ry) / (1 + ry) - c1 * (u[im+3*Mx+1] - u[im+1]) / (1 + ry);
		if(iz == Nl)
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[im+2] - u[2*Ndx+2*Ndy+3*(ix+iy*Mx)+3*(Nl-1)*Mxy+2]) / (1 + rz);
		else if(iz == Nz+Nl-1)
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[2*Ndx+2*Ndy+Ndz+3*(ix+iy*Mx)+2] - u[im-3*Nl*Mx+2]) / (1 + rz);
		else
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[im+2] - u[im-3*Nl*Mx+2]) / (1 + rz);
		dp[ip] = pm[im] + pm[im+1] + pm[im+2];
	}

// -z
	if(iz < Nl){
		im = 3 * (iz*Mxy + iy*Mx + ix) + 2*Ndx + 2*Ndy;
		rx = 0; ry = 0; rz = 0;
		if(ix <     Nl) rx = alpha * pow((float)(Nl-ix-0.5)/Nl, 2) / 2;
		if(ix >= Nx+Nl) rx = alpha * pow((float)(ix-Nx-Nl+0.5)/Nl, 2) / 2;
		if(iy <     Nl) ry = alpha * pow((float)(Nl-iy-0.5)/Nl, 2) / 2;
		if(iy >= Ny+Nl) ry = alpha * pow((float)(iy-Ny-Nl+0.5)/Nl, 2) / 2;
		rz = alpha * pow((float)(Nl-iz-0.5)/Nl, 2) / 2;
		if(ix > 0 && ix < Mx-1)
			pm[im]   = pm[im]   * (1 - rx) / (1 + rx) - c1 * (u[im] - u[im-3]) / (1 + rx);
		if(iy > 0 && iy < My-1)
			pm[im+1] = pm[im+1] * (1 - ry) / (1 + ry) - c1 * (u[im+1] - u[im-3*Mx+1]) / (1 + ry);
		if(iz > 0 && iz < Mz-1)
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[im+2] - u[im-3*Mxy+2]) / (1 + rz);
		dp[ip] = pm[im] + pm[im+1] + pm[im+2];
	}

// +z
	if(iz >= Nz+Nl && iz < Mz){
		im = 3 * ((iz-(Nz+Nl))*Mxy + iy*Mx + ix) + 2*Ndx + 2*Ndy + Ndz;
		rx = 0; ry = 0; rz = 0;
		if(ix <     Nl) rx = alpha * pow((float)(Nl-ix-0.5)/Nl, 2) / 2;
		if(ix >= Nx+Nl) rx = alpha * pow((float)(ix-Nx-Nl+0.5)/Nl, 2) / 2;
		if(iy <     Nl) ry = alpha * pow((float)(Nl-iy-0.5)/Nl, 2) / 2;
		if(iy >= Ny+Nl) ry = alpha * pow((float)(iy-Ny-Nl+0.5)/Nl, 2) / 2;
		rz = alpha * pow((float)(iz-(Nz+Nl)+0.5)/Nl, 2) / 2;
		if(ix > 0 && ix < Mx-1)
			pm[im]   = pm[im]   * (1 - rx) / (1 + rx) - c1 * (u[im] - u[im-3]) / (1 + rx);
		if(iy > 0 && iy < My-1)
			pm[im+1] = pm[im+1] * (1 - ry) / (1 + ry) - c1 * (u[im+1] - u[im-3*Mx+1]) / (1 + ry);
		if(iz < Mz-1)
			pm[im+2] = pm[im+2] * (1 - rz) / (1 + rz) - c1 * (u[im+3*Mxy+2] - u[im+2]) / (1 + rz);
		dp[ip] = pm[im] + pm[im+1] + pm[im+2];
	}

}


__global__ void PML_particle(float* dp, float* u, int3 Ndiv, int3 Ndim, int3 Nd, float c2, float alpha, int Nl){

	const unsigned int tx  = threadIdx.x;
	const unsigned int ty  = threadIdx.y;
	const unsigned int tz  = threadIdx.z;
	const unsigned int bdx = blockDim.x;
	const unsigned int bdy = blockDim.y;
	const unsigned int bdz = blockDim.z;
	const unsigned int bx  = blockIdx.x;
	const unsigned int by  = blockIdx.y;
	const unsigned int bz  = blockIdx.z;

	const unsigned int ix = bx * bdx + tx;
	const unsigned int iy = by * bdy + ty;
	const unsigned int iz = bz * bdz + tz;
	const unsigned int Nx = Ndim.x;
	const unsigned int Ny = Ndim.y;
	const unsigned int Nz = Ndim.z;
	const unsigned int Mx = Ndiv.x;
	const unsigned int My = Ndiv.y;
	const unsigned int Mz = Ndiv.z;
	const unsigned int Mxy = Ndiv.x*Ndiv.y;
	const unsigned int Ndx = 3*Nd.x;
	const unsigned int Ndy = 3*Nd.y;
	const unsigned int Ndz = 3*Nd.z;
	const unsigned int ip = iz * Mxy + iy * Mx + ix;

	long long im;
	float rx, ry, rz;


// -x
	if(ix < Nl && iy >= Nl && iy < Ny+Nl && iz >= Nl && iz < Nz+Nl){
		im = 3 * ((iz-Nl)*Nl*Ny + (iy-Nl)*Nl + ix);
		rx = 0; ry = 0; rz = 0;
		rx = alpha * pow((float)(Nl-ix-1)/Nl, 2) / 2;
		u[im] = u[im] * (1 - rx) / (1 + rx) - c2 * (dp[ip+1] - dp[ip]) / (1 + rx);
		if(iy < Ny+Nl-1) u[im+1] = u[im+1] * (1 - ry) / (1 + ry) - c2 * (dp[ip+Mx] - dp[ip]) / (1 + ry);
		if(iz < Nz+Nl-1) u[im+2] = u[im+2] * (1 - rz) / (1 + rz) - c2 * (dp[ip+Mxy] - dp[ip]) / (1 + rz);
	}

// +x
	if(ix >= Nx+Nl && ix < Mx && iy >= Nl && iy < Ny+Nl && iz >= Nl && iz < Nz+Nl){
		im = 3 * ((iz-Nl)*Nl*Ny + (iy-Nl)*Nl + ix-(Nx+Nl)) + Ndx;
		rx = 0; ry = 0; rz = 0;
		rx = alpha * pow((float)(ix-(Nx+Nl))/Nl, 2) / 2;
		u[im] = u[im]   * (1 - rx) / (1 + rx) - c2 * (dp[ip] - dp[ip-1]) / (1 + rx);
		if(iy < Ny+Nl-1) u[im+1] = u[im+1] * (1 - ry) / (1 + ry) - c2 * (dp[ip+Mx] - dp[ip]) / (1 + ry);
		if(iz < Nz+Nl-1) u[im+2] = u[im+2] * (1 - rz) / (1 + rz) - c2 * (dp[ip+Mxy] - dp[ip]) / (1 + rz);
	}

// -y
	if(iy < Nl && iz >= Nl && iz < Nz+Nl){
		im = 3 * ((iz-Nl)*Nl*Mx + iy*Mx + ix) + 2*Ndx;
		rx = 0; ry = 0; rz = 0;
		if(ix < Nl) rx = alpha * pow((float)(Nl-ix-1)/Nl, 2) / 2;
		if(ix >= Nx+Nl-1) rx = alpha * pow((float)(ix-Nx-Nl+1)/Nl, 2) / 2;
		ry = alpha * pow((float)(Nl-iy-1)/Nl, 2) / 2;
		if(ix < Mx-1) u[im] = u[im] * (1 - rx) / (1 + rx) - c2 * (dp[ip+1] - dp[ip]) / (1 + rx);
		u[im+1] = u[im+1] * (1 - ry) / (1 + ry) - c2 * (dp[ip+Mx] - dp[ip]) / (1 + ry);
		if(iz < Nz+Nl-1) u[im+2] = u[im+2] * (1 - rz) / (1 + rz) - c2 * (dp[ip+Mxy] - dp[ip]) / (1 + rz);
	}

// +y
	if(iy >= Ny+Nl && iz >= Nl && iz < Nz+Nl){
		im = 3 * ((iz-Nl)*Nl*Mx + (iy-(Ny+Nl))*Mx + ix) + 2*Ndx + Ndy;
		rx = 0; ry = 0; rz = 0;
		if(ix < Nl) rx = alpha * pow((float)(Nl-ix-1)/Nl, 2) / 2;
		if(ix >= Nx+Nl-1) rx = alpha * pow((float)(ix-Nx-Nl+1)/Nl, 2) / 2;
		ry = alpha * pow((float)(iy-(Ny+Nl))/Nl, 2) / 2;
		if(ix < Mx-1) u[im] = u[im] * (1 - rx) / (1 + rx) - c2 * (dp[ip+1] - dp[ip]) / (1 + rx);
		u[im+1] = u[im+1] * (1 - ry) / (1 + ry) - c2 * (dp[ip] - dp[ip-Mx]) / (1 + ry);
		if(iz < Nz+Nl-1) u[im+2] = u[im+2] * (1 - rz) / (1 + rz) - c2 * (dp[ip+Mxy] - dp[ip]) / (1 + rz);
	}

// -z
	if(iz < Nl){
		im = 3 * (iz * Mxy + iy*Mx + ix) + 2*Ndx + 2*Ndy;
		rx = 0; ry = 0; rz = 0;
		if(ix < Nl) rx = alpha * pow((float)(Nl-ix-1)/Nl, 2) / 2;
		if(ix >= Nx+Nl-1) rx = alpha * pow((float)(ix-Nx-Nl+1)/Nl, 2) / 2;
		if(iy < Nl) ry = alpha * pow((float)(Nl-iy-1)/Nl, 2) / 2;
		if(iy >= Ny+Nl-1) ry = alpha * pow((float)(iy-Ny-Nl+1)/Nl, 2) / 2;
		rz = alpha * pow((float)(Nl-iz-1)/Nl, 2) / 2;
		if(ix < Mx-1) u[im]   = u[im]   * (1 - rx) / (1 + rx) - c2 * (dp[ip+1] - dp[ip]) / (1 + rx);
		if(iy < My-1) u[im+1] = u[im+1] * (1 - ry) / (1 + ry) - c2 * (dp[ip+Mx] - dp[ip]) / (1 + ry);
		u[im+2] = u[im+2] * (1 - rz) / (1 + rz) - c2 * (dp[ip+Mxy] - dp[ip]) / (1 + rz);
	}

// +z
	if(iz >= Nz+Nl && iz < Mz){
		im = 3 * ((iz-(Nz+Nl))*Mxy + iy*Mx + ix) + 2*Ndx + 2*Ndy + Ndz;
		rx = 0; ry = 0; rz = 0;
		if(ix < Nl) rx = alpha * pow((float)(Nl-ix-1)/Nl, 2) / 2;
		if(ix >= Nx+Nl-1) rx = alpha * pow((float)(ix-Nx-Nl+1)/Nl, 2) / 2;
		if(iy < Nl) ry = alpha * pow((float)(Nl-iy-1)/Nl, 2) / 2;
		if(iy >= Ny+Nl-1) ry = alpha * pow((float)(iy-Ny-Nl+1)/Nl, 2) / 2;
		rz = alpha * pow((float)(iz-(Nz+Nl))/Nl, 2) / 2;
		if(ix < Mx-1) u[im]   = u[im]   * (1 - rx) / (1 + rx) - c2 * (dp[ip+1] - dp[ip]) / (1 + rx);
		if(iy < My-1) u[im+1] = u[im+1] * (1 - ry) / (1 + ry) - c2 * (dp[ip+Mx] - dp[ip]) / (1 + ry);
		u[im+2] = u[im+2] * (1 - rz) / (1 + rz) - c2 * (dp[ip] - dp[ip-Mxy]) / (1 + rz);
	}

}


// 音圧分布保存
void save_cross_section(float* dp, float* pp, int it, char* DirName)
{
	float* pobs  = (float*) malloc(sizeof(float));		 	// 観測点音圧
	int id, idd;
	char filename[200];
	const unsigned int Mx = Ndiv.x;
	const unsigned int My = Ndiv.y;
	const unsigned int Mz = Ndiv.z;

	// xy平面保存
	if(iplane == 1){
		for(int j = 0; j < My; j++){
			for(int i = 0; i < Mx; i++){
				id = ipn * Mx * My + j * Mx + i;
				idd = j * Mx + i;
				cudaMemcpy(pobs, dp+id, sizeof(float), cudaMemcpyDeviceToHost);
				pp[idd] = (*pobs);
			}
		}
		if(iwave == 0){
			sprintf(filename, "%soutput%d.dat", DirName, it);	// ボクセルファイル
			FILE *fpo  = fopen(filename,"w");
			for(int j = 0; j < My; j++){
				for(int i = 0; i < Mx; i++){
					id = j * Mx + i;
					fprintf(fpo, "%e\t", pp[id]);
				}
				fprintf(fpo, "\n");
			}
			fclose(fpo);
		}
		else{
			sprintf(filename, "%soutput%d.bin", DirName, it);	// 出力ファイル
			FILE *fpo  = fopen(filename,"wb");
			fwrite(&Mx, sizeof(int), 1, fpo);
			fwrite(&My, sizeof(int), 1, fpo);
			fwrite(pp, sizeof(float), Mx*My, fpo);
			fclose(fpo);
		}
	}

	// yz平面保存
	if(iplane == 2){
		for(int k = 0; k < Mz; k++){
			for(int j = 0; j < My; j++){
				id = k * Mx * My + j * Mx + ipn;
				idd = k * My + j;
				cudaMemcpy(pobs, dp+id, sizeof(float), cudaMemcpyDeviceToHost);
				pp[idd] = (*pobs);
			}
		}
		if(iwave == 0){
			sprintf(filename, "%soutput%d.dat", DirName, it);	// ボクセルファイル
			FILE *fpo  = fopen(filename,"w");
			for(int k = 0; k < Mz; k++){
				for(int j = 0; j < My; j++){
					id = k * My + j;
					fprintf(fpo, "%e\t", pp[id]);
				}
				fprintf(fpo, "\n");
			}
			fclose(fpo);
		}
		else{
			sprintf(filename, "%soutput%d.bin", DirName, it);	// 出力ファイル
			FILE *fpo  = fopen(filename,"wb");
			fwrite(&My, sizeof(int), 1, fpo);
			fwrite(&Mz, sizeof(int), 1, fpo);
			fwrite(pp, sizeof(float), My*Mz, fpo);
			fclose(fpo);
		}
	}

	// xz平面保存
	if(iplane == 3){
		for(int k = 0; k < Mz; k++){
			for(int i = 0; i < Mx; i++){
				id = k * Mx * My + ipn * Mx + i;
				idd = k * Mx + i;
				cudaMemcpy(pobs, dp+id, sizeof(float), cudaMemcpyDeviceToHost);
				pp[idd] = (*pobs);
			}
		}
		if(iwave == 0){
			sprintf(filename, "%soutput%d.dat", DirName, it);	// ボクセルファイル
			FILE *fpo  = fopen(filename,"w");
			for(int k = 0; k < Mz; k++){
				for(int i = 0; i < Mx; i++){
					id = k * Mx + i;
					fprintf(fpo, "%e\t", pp[id]);
				}
				fprintf(fpo, "\n");
			}
			fclose(fpo);
		}
		else{
			sprintf(filename, "%soutput%d.bin", DirName, it);	// 出力ファイル
			FILE *fpo  = fopen(filename,"wb");
			fwrite(&Mx, sizeof(int), 1, fpo);
			fwrite(&Mz, sizeof(int), 1, fpo);
			fwrite(pp, sizeof(float), Mx*Mz, fpo);
			fclose(fpo);
		}
	}

}


// 配列表示
void print_matrix(unsigned char* Vox, int ipx1, int ipx2, int ipy1, int ipy2, int ipz1, int ipz2)
{
	int i, j, k;
	unsigned long long id;
	
	if(ipx2 > Ndiv.x) ipx2 = Ndiv.x;
	if(ipy2 > Ndiv.y) ipy2 = Ndiv.y;
//	if(ipz2 > Ndiv.z) ipz2 = Ndiv.z;
	for(k = ipz2-1; k >= ipz1; k--){
		printf("k=%d\n",k);
		for(j = ipy2-1; j >= ipy1; j--){
			for(i = ipx1; i < ipx2; i++){
				id = Ndiv.x * Ndiv.y * k + Ndiv.x * j + i;
				printf("%2d ", Vox[id]);
			}
			printf("\n");
		}
		printf("\n");
	}
}


__global__ void DivFinder(float* dp, int Nx, int Ny, int Nzdiv, float* dMax, 
	unsigned long long* dMid, int iReg)
{
	
	int ix, iy, iz;
	float ap;
	unsigned long long id;

	dMax[iReg] = 0;
	dMid[iReg] = 0;
	for(iz = 0; iz < Nzdiv; iz++){
		for(iy = 0; iy < Ny; iy++){
			for(ix = 0; ix < Nx; ix++){
				id = Nx * Ny * iz + Nx * iy + ix;
				ap = abs(dp[id]);
				if(ap > dMax[iReg]){
					dMax[iReg] = ap;
					dMid[iReg] = id;
				}
			}
		}
	}
}


void DivError(float* dp, int Nzdiv, int Nreg, int iReg)
{
	double* hMax = (double*) malloc(sizeof(double)*Nreg);		// 観測点発散値
	unsigned long long* hMid = (unsigned long long*) malloc(sizeof(unsigned long long)*Nreg);	// 観測点発散id
	float *dMax;
	unsigned long long *dMid;

	for(int i = 0; i < Nreg; i++){
		hMax[i] = 0;
		hMid[i] = 0;
	}
	cudaMalloc((void**) &dMax, sizeof(float)*Nreg);					// 発散値
	cudaMalloc((void**) &dMid, sizeof(unsigned long long)*Nreg);	// 発散id
	cudaMemset(dMax,  0, sizeof(float)*Nreg);
	cudaMemset(dMid,  0, sizeof(unsigned long long)*Nreg);

	DivFinder<<<1, 1>>>(dp, Ndiv.x, Ndiv.y, Nzdiv, dMax, dMid, iReg);
	cudaMemcpy(hMax, dMax, sizeof(float)*Nreg, cudaMemcpyDeviceToHost);
	cudaMemcpy(hMid, dMid, sizeof(unsigned long long)*Nreg, cudaMemcpyDeviceToHost);
	for(int ir = 0; ir < Nreg; ir++)
		printf("iReg: %d id=%lld: %e\n", ir, hMid[ir], hMax[ir]);
}
