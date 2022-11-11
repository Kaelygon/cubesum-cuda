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

//KERNEL FUNCTIONS
//(for any natural a,b a>b tdtarg3=a^3+b^3+c^3 find c by quadratic convergence)
__global__
void nmcbrt(__uint64_t n, __uint128_t tdtarg3, __uint64_t tdtarg, __uint128_t aofs, vdint3 **result ){
	__uint64_t index = blockIdx.x * blockDim.x + threadIdx.x;

	__uint64_t stride = blockDim.x * gridDim.x;

	__uint64_t a,b,c;
	__uint128_t ai3,bc3,ctarg;
	__uint128_t c0=tdtarg-1;
	//search a,b. Index starts from 1
	for (a = index+aofs+1; a<n+aofs+1; a += stride){
		int inc=0;
		ai3=(__uint128_t)a*a*a;
		bc3=tdtarg3-ai3; //b^3+c^3

		c0=tdtarg-1;

		//quadratic iterations c for (a^3+(a-1)^3+c^3)^(1/3)
		ctarg=bc3-ai3;
		do{ 
			c = c0;
			c0 = ((__uint128_t)3*c + ctarg/((__uint128_t)c*c))/4; //newton's method
//			c0 = ((3*ctarg)/((2*c*c*c+ctarg))-1)/(2*c)+c*c*c-ctarg; //halley's steps
		}while (c0 < c);
		
		for (b = a+1; b<c; b++) {
			ctarg=bc3-(__uint128_t)b*b*b; //c^3

			//only one newton's method pass is required for each b++.
			c = ((__uint128_t)3*c + ctarg/((__uint128_t)c*c))/4;

			//if ctarg^(1/3) has integer soltion
			((__uint128_t)c*c*c==ctarg) ? result[inc++][a-aofs-1]=vdint3{a,b,c} : vdint3() ;
			#if DEBUG_DIGS==1
				vdint3{a,b,c}.printvec();
			#endif
		}
	}
}

//use builtin cbrt 
__global__
void bicbrt(__uint64_t n, __uint128_t tdtarg3, __uint64_t tdtarg, __uint128_t aofs, vdint3 **result ){
	__uint64_t index = blockIdx.x * blockDim.x + threadIdx.x;

	__uint64_t stride = blockDim.x * gridDim.x;

	__uint64_t a,b,c;
	__uint128_t bc3,ctarg;
	
	//search a,b. Index starts from 1
	for (a = index+aofs+1; a<n+aofs+1; a += stride){
		int inc=0;
		bc3=tdtarg3-(__uint128_t)a*a*a; //b^3+c^3

		c=tdtarg-1;
		for (b = a+1; b<c; b++) {
			ctarg=bc3-(__uint128_t)b*b*b; //c^3

			c=__builtin_cbrt((double)ctarg);

			//if ctarg^(1/3) has integer soltion
			((__uint128_t)c*c*c==ctarg) ? result[inc++][a-aofs-1]=vdint3{a,b,c} : vdint3() ;

			#if DEBUG_DIGS==1
				vdint3{a,b,c}.printvec();
			#endif
		}
	}
}
//search every number of sequence A023042 up to tdtarg
__global__
void allcbrt(__uint64_t n, __uint128_t tdtarg3, __uint64_t tdtarg, __uint128_t aofs, vdint3 **result ){
	__uint64_t index = blockIdx.x * blockDim.x + threadIdx.x;

	__uint64_t stride = blockDim.x * gridDim.x;

	__uint64_t a,b,c,cbi;
	__uint128_t ai3,bi3,bc3,ctarg;
	__uint128_t c0=tdtarg-1;

	__uint128_t cube;
	for(cbi=index+aofs;cbi<=n+aofs;cbi+=stride){ //in this case, cbi is the target  
		cube=cbi*cbi*cbi;
		for (a = 1; a<cbi*0.69336+1; a += 1){
			int inc=0;
			ai3=(__uint128_t)a*a*a;
			if(cube<ai3){break;}
			bc3=cube-ai3; //b^3+c^3

			c0=tdtarg-1;

			//quadratic iterations c for (a^3+(a-1)^3+c^3)^(1/3)
			ctarg=bc3-ai3;
			do{ 
				c = c0;
				c0 = ((__uint128_t)3*c + ctarg/((__uint128_t)c*c))/4; //newton's method
			}while (c0 < c);
			
			for (b = a+1; b<c; b++) {
				bi3=b*b*b;
				if(bc3<bi3){break;}
				ctarg=bc3-(__uint128_t)bi3; //c^3

				//only one newton's method pass is required for each b++.
				c = ((__uint128_t)3*c + ctarg/((__uint128_t)c*c))/4;

				#if DEBUG_DIGS==0
					vdint3{a,b,c}.printvec();
				#endif

				//if ctarg^(1/3) has integer soltion
				if((__uint128_t)c*c*c==ctarg){
					result[inc++][a]=vdint3{a,b,c};
					goto iend;
				}
			}
		}
		iend:;
	};
}

//array of algorithm functions. It's like the entire internet hasn't heard about these
static void  (*algoarr[3])(__uint64_t n, __uint128_t tdtarg3, __uint64_t tdtarg, __uint128_t aofs, vdint3 **result ) 
	= {nmcbrt, bicbrt, allcbrt};

//EOF KERNEL FUNCTIONS

void searchResults(const int tdid, const int tdcount, const int rpt, __uint128_t tdtasks, const string tfil, vdint3 **result, char *resultstr, int *tdfound){
	string tmpstr="";
	for( __uint128_t k=tdid; k<tdtasks; k+=tdcount ){
		if(result[0][k].x==0){continue;}
		for(__uint128_t j=0;j<rpt;j++ ){
			if(result[j][k].x==0){continue;}
			tmpstr+= result[j][k].stringvec()+"\n";
			*tdfound=*tdfound+1;
			result[j][k]=vdint3(); //reset value
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
	getline(cin, input);
	if ( input=="q" || input=="quit" )
	{
		std::cout << "#Forced to quit\n";
		exit(-1);
	}
	return;
}

int main(int argc, char **argv){
	/*CLOCK*/double d_time;
	
	//BOF config
	
	uint tc=1; //thread count
	__uint64_t targ = 131071;
	__uint64_t update_rate = targ; //gpu workload size before going back to host, tasksize=targ for best performance, but no logging of progress
	uint maxrpb=2; //how many results with same A value can be stored. I have only found values which have two solutions with same A, e.g. targ= 189(odd), 256(even), ... 6959(prime)
	__uint64_t start = 0; //A offset

	string work_directory = "./";
	string results_file = "results.txt";
	string config = work_directory+"config.cfg";	//config file
	string progress_file = work_directory+"lastSolve.txt";
	__uint64_t maxvram=1*1024*1024*1024;
	uint algo=0;

	bool clear_file=0;
	
	if( argc>1 ){
		config = argv[1];
	}

	ifstream config_file (config);
	if(config_file.good()){
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
					if(start!=0){start-=1;} //search starts from start+1
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
				}else if(	word=="algorithm" ){
					algo=stol(value);
				}
				cout << "#" << word << "=" << value << "\n";
				word="";
			}
		}
	}else{
		printf("Config not found.\n");
		exit(-1);
	}
	config_file.close();

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

	__uint128_t targ3 = (__uint128_t)targ*targ*targ;

	__uint64_t tasks = 
	(algo==2) ? 
		targ-start 	
		:
		(targ)*(__float128)0.693361274350634659846548402128973976+1-start //max A = t * 3^(2/3)/3
	;

	__uint64_t taskblocks;
	if(tc>tasks){tc=tasks;} //limit max threads to tasks

	//calculate vram. if max is exceeded: recalc
	__uint64_t arrx = maxrpb*sizeof(vdint3*);
	__uint64_t arry = maxrpb*tasks*sizeof(vdint3);	// algo==2 has same allocation size as there's cbrt(targ) perfect cubes up to targ
	
	__uint64_t totalloc = arrx + arry;
	if( totalloc > maxvram ){
		taskblocks=totalloc/maxvram+1;
		tasks=(tasks-1)/taskblocks+1;
		arry = maxrpb*tasks*sizeof(vdint3);
		totalloc = arrx + arry;
		cout << "New update rate due to vram limit: " << totalloc << " taskblocks: " << taskblocks << "\n";
	}else{
		taskblocks=(tasks-1)/update_rate+1;
		tasks=(tasks-1)/taskblocks+1;
	}
	cout << "#Allocated: " << totalloc/1024 << " KiB\n";

	vdint3 **result; //shared with host and device

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
	__uint64_t lastsolve=start;

	ofstream out_result_file;
	out_result_file.open(work_directory+results_file, ios::out | ios::app);

	//splitting up gpu progress
	for(__uint64_t ti=0;ti<taskblocks;ti++){

		__uint64_t blockSize = 1024;
		if(blockSize>tasks){blockSize=tasks/32*32+32;}	//prevent initializing beyond tasks, int truncate 32s
		__uint64_t numBlocks = ((tasks-1)/taskblocks) / blockSize + 1;

		//GPU
		algoarr[algo]<<<numBlocks,blockSize>>>(
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
			
			out_result_file << (char*)resultstr[tid];
			cout << (char*)resultstr[tid];

			free(resultstr[tid]);
		}
		
		delete [] wres;
		free(tdfound);
		free(resultstr);
		resultstr=NULL;
		//EOF wres

		__uint64_t curprog = start+tasks*(ti+1)+1; //current (progress of) A
		out_result_file << "#a: "+ui128tos( curprog )+"\n" << flush;
		cout << "#a: "+ui128tos( curprog )+"\n" << flush;

		lastsolve=curprog;
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
//	cudaDeviceReset();

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
	
	{//approx iter/s scope
		__int128_t aprx = targ*0.69336-start-(targ*0.69336-lastsolve-2); // [checked A count] - start - ([not checked A count]) 
		long double iters = aprx*aprx /4.5 /(long double)d_time;
		string magsfx[5] = {""," Kilo"," Million"," Billion"," Giga"};
		int mi=0;
		while(iters>1000 && mi<5){
			iters/=1000;
			mi++;
		}
		if((algo==0 || algo==1)){
			/*CLOCK*/cout << "#" << iters << magsfx[mi] << " C per second\n";
		}
	}
	
	return 0;
}
// nvprof ./cubesum
// /opt/cuda/extras/compute-sanitizer/compute-sanitizer --leak-check full ./cubesum
// cuda-memcheck ./cubesum