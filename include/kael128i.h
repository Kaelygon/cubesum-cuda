#define MAX__UINT128 (__uint128_t)(-1)
#define MAX__UINT64 (__uint64_t)(-1)

using namespace std;

//128 - (count leading zeros ui128) + 2 
inline int iclz_ui128 (__uint128_t u) {
	uint64_t hi = u>>64;
	uint64_t lo = u;

	return  
		hi ? 
		 	lo&hi ?  66-__builtin_clzll(lo) : 2								//hi true  -> lo&hi ? 1 : 0
			:
			lo&hi ? 130-__builtin_clzll(hi) : 66-__builtin_clzll(lo)		//hi false -> lo&hi ? 2 : 1
	;
}

//unsigned integer quadruple cube root
__uint128_t uiqcbrt(__uint128_t n){

	__uint128_t r;
	int lz=iclz_ui128(n)/3;
	__uint128_t r0=1<<lz;

	do{ //quadratic iterations
		r = r0;
		r0 = (3*r + n/(r*r))/4 ;
	}
	while (r0 < r);

	return (r*r*r==n) ? r : 0;
		
}

template <class T>
string ui128tos(T n){//128 uint to string
	if(n==0){return (string)"0";}
	
	string str = "";

	while(n){
		str+=to_string( (__uint8_t)(n%10) );
		n/=10;
	}
	
	string inv = "";
	for( int i=str.size()-1;i>=0;i-- ){
		inv+=str[i];
	}

	return inv;
}

//unsigned quadruple (precision) integer [3]
struct uqint3
{
    __uint128_t x, y, z;
};
struct vqint3
{
    __uint128_t x, y, z;

    __host__ __device__ constexpr vqint3(__uint128_t qx = 0, __uint128_t qy = 0, __uint128_t qz = 0) : x(qx), y(qy), z(qz) {}
    __host__ __device__ constexpr vqint3(uqint3 q) : x(q.x), y(q.y), z(q.z) {}
    __host__ __device__ inline void printvec() {
		printf(
			"%lld^3 + %lld^3 + %lld^3\n",
			(unsigned long long)(x & 0xFFFFFFFFFFFFFFFF),
			(unsigned long long)(y & 0xFFFFFFFFFFFFFFFF),
			(unsigned long long)(z & 0xFFFFFFFFFFFFFFFF)
		);
	}
    __host__ string stringvec() { //max 134 char
		return
			ui128tos(x) + "^3" + " + " +
			ui128tos(y) + "^3" + " + " +
			ui128tos(z) + "^3" 
		;
	}

};