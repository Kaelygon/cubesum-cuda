Nvidia Cuda variant of the program "cubesum" at https://github.com/Kaelygon/cubesum

Finds N with property that N^3=A^3+B^3+C^3 multithreaded 
OEIS A023042 https://oeis.org/A023042/
C>B>A>=0 N>2

Usage ./cubesum [config file]

If no argument is given, default file is ./config.cfg, if file or some elements dont exist the program uses following default values:

	thread_count=24			# (uint) thread count
	start=0				# (__uint64_t)
	target=262143			# (__uint64_t) target number N.
	update_rate=0			# (__uint64_t) split gpu work into [update_rate] sized chunks 0=disabled
	work_directory=./datacs/	# (string)
	results_file=results.txt	# (string) matches stored here
	progress_file=#lastSolve.txt	# (string) progress file, if any
	clear_file=true			# (bool) clear files at start?
	max_vram=6G			# (__uint64_t) max allocated video memory. Suffixes ""=Bytes, "K"=KiB, "M"=MiB, "G"=GiB 
	results_per_block=2		# (uint) maximum stored results with sane A per block
	algorithm=1			# (uint) 0 = (Newton's method) 1 = (__builtin_cbrt)
	?				# end of file

Elements order don't matter in config file. Comment lines with '#'
Question mark (?) sets end of file

You can enter q or quit during run and the program will save last solved up to n * A * update_rate, to [work_directory]/[progress_file]
Enter q again to force quit without computing current split

Macros to compile

Use Newton's method
 ALGO 0
 max target = (2^128)^(1/3) ~= 2^46.7

Use builtin cbrt (64-bit double)
 ALGO 1
  Fast but precision is limited to (2^52-1)^(1/3) = 208063 
  C higher than this may not compute properly in some special cases
  Cuda/gcc does not support float128 cbrt as of now.