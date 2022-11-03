#include <ctime>
	/*CLOCK*/clock_t st_time = clock();

#include <iostream>
#include <thread>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

#include "./include/kael128i.h"

//https://github.com/NVIDIA/cuda-samples/tree/master/Common
#define DEBUG_KAEL 0
#if DEBUG_KAEL==1
	#include <helper_cuda.h> 
#endif

// Kernel function, (for any natural a,b a>b tdtarg3=a^3+b^3+c^3 find c by quadratic convergence)
__global__
void quadCverg(__uint128_t n, __uint128_t tdtarg3, __uint128_t tdtarg, __uint128_t aofs, vqint3 **result ){
	__uint128_t index = blockIdx.x * blockDim.x + threadIdx.x;

	__uint128_t stride = blockDim.x * gridDim.x;

	__uint128_t a,b,c;
	__uint128_t ai3,ctarg;
	__uint128_t c0=tdtarg-1; //c0 value is reused in each loop of b
	
	//search a,b. Index starts from 1
	for (a = index+aofs+1; a<n+aofs; a += stride){
		int inc=0;
		ai3=a*a*a;

		c=tdtarg-1;
		c0=c;
		for (b = a+1; b<c; b++) {
			ctarg=tdtarg3-ai3-b*b*b; //ctarg==c*c*c

			do{ //quadratic iterations c = ctarg-(c^3-ctarg)/derivative(c^3-ctarg,c)
				c = c0;
				c0 = (3*c + ctarg/(c*c))/4;
			}while (c0 < c);
				
			//save or print result
			(c*c*c==ctarg) ? result[inc++][a-aofs-1]=uqint3{a,b,c} : vqint3() ;

			#if DEBUG_DIGS==1
				vqint3{a,b,c}.printvec();
			#endif
		}
	}
}

void searchResults(const int tdid, const int tdcount, const int rpt, __uint128_t n, const string tfil, vqint3 **result, char *resultstr, int *tdfound){
	string tmpstr="";
	for( __uint128_t k=tdid; k<n; k+=tdcount ){
		if(result[0][k].x==0){continue;}
		for(__uint128_t j=0;j<rpt;j++ ){
			if(result[j][k].x==0){continue;}
			tmpstr+= result[j][k].stringvec()+"\n";
			*tdfound=*tdfound+1;
			result[j][k]=vqint3(); //reset value
		}
	}
	strcpy(resultstr,tmpstr.c_str());
	return;
}

bool threads_active=1;
void bgtasks(void)
{
	string input;
	while(threads_active==1){
		getline(cin, input);
		if ( input=="q" || input=="quit" )
		{
			std::cout << "#User exit\n";
			std::cout << "#Closing threads\n";
			threads_active=0;
		}
	}
	return;
}

int main(int argc, char **argv){
	/*CLOCK*/double d_time;
	
	//BOF config
	
	int tc=1; //thread count
	__uint128_t targ = 131071;
	__uint128_t update_rate = targ; //gpu workload size before going back to host, tasksize=targ for best performance, but no logging of progress
	int maxrpb=2; //how many results with same A value can be stored. I have only found values which have two solutions with same A, e.g. targ= 189(odd), 256(even), ... 6959(prime)
	__uint128_t start = 0; //A offset

	string work_directory = "./";
	string results_file = "results.txt";
	string config = work_directory+"config.cfg";	//config file
	string progress_file = work_directory+"lastSolve.txt";
	__uint128_t maxvram=1*1024*1024*1024;
	bool clear_file=0;
	
	if( argc>1 ){
		config = argv[1];
	}

	if(config!=""){
		ifstream config_file (config);
		stringstream tmp;
		tmp << config_file.rdbuf();
		string config_string = tmp.str();

		//parser
		string word="";
		for(uint i=0; (i<config_string.size()) && (config_string[i]!='?') ;i++){ //word before '='
			if(config_string[i]!='='){
				if(config_string[i]==' ' || config_string[i]=='	'){continue;} //ignore white spaces
				word+=config_string[i];
				if(config_string[i]=='#'){//skip comments
					while(config_string[i]!='\n' && i<config_string.size() ){i++;}word="";
				}
			}else{
				string value="";
				i++;
				for(i=i;i<config_string.size();i++){ //value after '=' or '#'
					if(config_string[i]==' ' || config_string[i]=='	'){continue;} //ignore white spaces
					if(config_string[i]=='#'){//skip comments
						while(config_string[i]!='\n'){i++;}break;
					}
					if(config_string[i]!='\n'){
						value+=(char)config_string[i];
					}else{break;}
				}

				//only few elements so this unelegant solution will do. std::map?
				if		(	word=="thread_count" ){
					tc=stol(value); 
					if(tc==0){tc=1;}
				}else if(	word=="start" ){
					start=stol(value); 
				}else if(	word=="target" ){
					targ=stol(value); 
					if(targ<=0){targ=1;}
				}else if(	word=="progress_file"){
					progress_file=value;
				}else if(	word=="results_file" ){
					results_file=value;
				}else if(	word=="update_rate" ){
					update_rate=stol(value); 
					if(update_rate==0){update_rate=targ;}
				}else if(	word=="work_directory" ){
					work_directory=value;
				}else if(	word=="clear_file" ){
					if(value=="true"){
						clear_file=1;	
					}else if(value=="false"){
						clear_file=0;	
					}else{
						clear_file=stol(value);
					}
				}else if(	word=="results_per_block" ){
					maxrpb=stol(value);
				}else if(	word=="max_vram" ){
					char vram_sfx=value[value.size()-1];
					if(vram_sfx == 'G' || vram_sfx == 'g' ){
						value.pop_back();
						maxvram=stol(value)*pow(1024,3);
						value+='G';
					}else if(vram_sfx == 'M' || vram_sfx == 'm' ){
						value.pop_back();
						maxvram=stol(value)*pow(1024,2);
						value+='M';
					}else if(vram_sfx == 'K' || vram_sfx == 'k' ){
						value.pop_back();
						maxvram=stol(value)*pow(1024,1);
						value+='K';
					}else{
						maxvram=stol(value);
					}
				}
				cout << "#" << word << "=" << value << "\n";
				word="";
			}
		}
	}

	if(clear_file==1){
		ofstream resfil;
		resfil.open(work_directory+results_file, ios::out);
		resfil << "";
		resfil.close();

		ofstream out_progress_file;
		out_progress_file.open(work_directory+progress_file, ios::out);
		out_progress_file << "0";
		out_progress_file.close();
	}
	
	if(progress_file!=""){//progress file
		ifstream in_progress_file;
		in_progress_file.open(work_directory+progress_file, ios::in);

		stringstream tmp;
		tmp << in_progress_file.rdbuf();
		in_progress_file.close();

		string progress_string = tmp.str();

		string word="";
		for(int i=0;progress_string[i];i++){ //word before '='
			if(progress_string[i]!='\n'){
				word+=progress_string[i];
			}
		}
		if(word!=""){
			start=stol(word);
		}
		cout << "#progress_file=" << to_string((uint64_t)start) << "\n";
	}

	//EOF config

	thread bgthread;
	bgthread = thread(bgtasks);

	__uint128_t targ3 = targ*targ*targ;

	__uint128_t tasks=targ*(__float128)0.693361274350634659846548402128973976+1; //max A = t * 3^(2/3)/3
	__uint128_t taskblocks;
	if(tc>tasks){tc=tasks;} //limit max threads to tasks

	//calculate vram. if max is exceeded: recalc
	__uint64_t arrx = maxrpb*sizeof(vqint3*);
	__uint64_t arry = maxrpb*tasks*sizeof(vqint3);
	__uint64_t totalloc = arrx + arry;
	if( totalloc > maxvram ){
		taskblocks=totalloc/maxvram+1;
		tasks=(tasks-1)/taskblocks+1;
		arry = maxrpb*tasks*sizeof(vqint3);
		totalloc = arrx + arry;
		cout << "New update rate due to vram limit: " << totalloc << "\n";
	}else{
		taskblocks=(tasks-1)/update_rate+1;
		tasks=(tasks-1)/taskblocks+1;
	}
	cout << "#Allocated: " << totalloc/1024 << " KiB\n";

	vqint3 **result; //shared with host and device

	#if DEBUG_KAEL==1
		checkCudaErrors( cudaMallocManaged( &result, arrx ) );
		for(int i=0;i<maxrpb;i++){
			checkCudaErrors( cudaMallocManaged( &result[i], arry ) );
		}
	#else
		cudaMallocManaged( &result, arrx );
		for(int i=0;i<maxrpb;i++){
			cudaMallocManaged( &result[i], arry );
		}
	#endif

	__uint64_t foundsum=0;
	__uint128_t lastsolve=start;

	ofstream out_result_file;
	out_result_file.open(work_directory+results_file, ios::out | ios::app);

	//splitting up gpu progress
	for(__uint128_t ti=start/tasks;ti<taskblocks;ti++){

		__uint128_t blockSize = 1024;
		if(blockSize>tasks){blockSize=tasks/32*32+32;}	//prevent initializing beyond tasks, int truncate 32s
		__uint128_t numBlocks = (tasks-1) / blockSize + 1;

		//GPU
		quadCverg<<<numBlocks,blockSize>>>(
			tasks, 
			targ3, 
			targ, 
			start+tasks*ti,
			result
		);
		
		// Wait for GPU to finish before accessing on host
		#if DEBUG_KAEL==1
			checkCudaErrors(cudaDeviceSynchronize());
		#else
			cudaDeviceSynchronize();
		#endif


		//THREADS search and write gpu results. wres
		int *tdfound = (int*) calloc(tc, sizeof(int*));
		char **resultstr = (char**) malloc(tc * sizeof(char*));		

		thread *wres = new thread[tc];

		for(uint tid=0;tid<tc;tid++ ){
			resultstr[tid] = (char*) malloc( (maxrpb*tasks)/tc * sizeof(char)*(135+5));
			wres[tid] = thread(
				searchResults,
					tid,
					tc,
					maxrpb,
					tasks,
					work_directory+results_file,
					&*result,
					&*resultstr[tid],
					&tdfound[tid]
			);
		}

		for(int tid=0;tid<tc;tid++){
			//synchronize threads
			wres[tid].join();
			foundsum+=tdfound[tid];
			
			out_result_file << resultstr[tid];

			free(resultstr[tid]);
		}
		
		delete [] wres;
		free(tdfound);
		free(resultstr);
		resultstr=NULL;
		//EOF wres

		out_result_file << "a: "+ui128tos( tasks*(ti+1) )+"\n";

		lastsolve=tasks*(ti+1);
		if(threads_active==0){
			break;
		}
	}
	out_result_file.close();
	
	#if DEBUG_KAEL==1
		for(int i=0;i<maxrpb;i++){
			checkCudaErrors( cudaFree(result[i]) );
		}
		checkCudaErrors( cudaFree(result) );
	#else
		for(int i=0;i<maxrpb;i++){
			cudaFree(result[i]);
		}
		cudaFree(result);	
	#endif
	cudaDeviceReset();

	ofstream out_progress_file;
	out_progress_file.open(work_directory+progress_file, ios::out);
	out_progress_file << ui128tos(lastsolve);
	out_progress_file.close();

	printf("\n#%s\n", cudaGetErrorString(cudaGetLastError()));
	printf("#Last solved: %lld\n",(__uint64_t)lastsolve);
	printf("#Found: %d\n",foundsum);


	bgthread.detach(); //will not free memory

	/*CLOCK*/d_time = (float)(clock()-st_time)/1000000;
	/*CLOCK*/cout << "#" << d_time << " s\n";
	/*CLOCK*/cout << "#" << (targ-start)*(targ-start)/7/1000000000/(long double)d_time << " Bil iter/s\n"; //appoximation

	return 0;
}
// nvprof ./cubesum
// /opt/cuda/extras/compute-sanitizer/compute-sanitizer --leak-check full ./cubesum
// cuda-memcheck ./cubesum