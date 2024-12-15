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
3
//Migration matrix 0
0 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
//Migration matrix 1
0 0 0 mig03
0 0 0 0
0 0 0 0
mig30 0 0 0
//Migration matrix 2
0 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
7 historical event
Tmigstop 0 0 0 1 0 0
Tmigstart 0 0 0 RESIZE4 0 1
TAdm 2 2 1 RESIZE2 0 2
TAdm 2 1 a_1_2 1 0 2
TAdm 2 0 1 1 0 2
Tdiv2 3 1 1 RESIZE3 0 2
Tdiv1 0 1 1 RESIZE1 0 2
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 1e-7 OUTEXP
