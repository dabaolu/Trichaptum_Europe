//Parameters for the coalescence simulation program : simcoal.exe
4 samples to simulate :
//Population effective sizes (number of genes)
Npop0
Npop1
Npop2
Npop3
//Samples sizes and samples age 
13
13
13
13
//Growth rates: negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
7 historical event
TAdm2 3 3 1 RESIZE3 0 0
TAdm2 3 1 a_1_3 1 0 0
TAdm2 3 0 1 1 0 0
TAdm1 2 2 1 RESIZE2 0 0
TAdm1 2 1 a_1_2 1 0 0
TAdm1 2 0 1 1 0 0
TDIV01 0 1 1 RESIZE1 0 0
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 1e-7 OUTEXP
