Nvidia Cuda variant of the program "cubesum" at https://github.com/Kaelygon/cubesum

Finds N with property that N^3=A^3+B^3+C^3 multithreaded 
OEIS A023042 https://oeis.org/A023042/
C>B>A>=0 N>2

Usage ./cubesum [config file]

If no argument is given, default file is ./config.cfg, if file or some elements dont exist the program uses following default values:

	thread_count=24			# (int) thread count
	start=0				# (__uint64_t)
	target=262143			# (__uint64_t) target number N.
	results_file=results.txt	# (string) matches stored here
	progress_file=#lastSolve.txt	# (string) progress file, if any
	update_rate=0			# (__uint64_t) split gpu work into [update_rate] sized chunks 0=disabled
	work_directory=./datacs/	# (string)
	clear_file=true			# (bool) clear files at start?
	max_vram=6G			# (__uint64_t) max allocated video memory
	results_per_block=2		# (int) maximum stored A results pre block
	?				# end of file

Elements order don't matter in config file. Comment lines with '#'
Question mark (?) sets end of file

You can enter q or quit during run and the program will save last solved up to n * A * update_rate, to [work_directory]/[progress_file]
Enter q again to force quit without computing current split

Some values are __uint128_t as target^3 exceeds ULL at 2^46.7

Macros to compile

Use builtin cbrt 
 USE_BI_CBRT 0
 max target = (2^128)^(1/3) ~= 2^46.7

 USE_BI_CBRT 1
  Fast but precision is limited to (2^52-1)^(1/3) = 208063 
  Meaning values with higher C than this with no trailing zeros, might not compute properly as 208063^3>2^52
  Cuda/gcc does not support float128 cbrt as of now.