typedef struct{										//creates struct that keeps track of
	char element[3];									//an element's data including:
	int val;											//name, number of valence electrons, and the electronegativity
	float e_neg;
	int a_num;
	float dist;
	float radius;
} elem_data;

elem_data table[16] = {							//array that keeps track of the name, num of val electrons and electronegativity
	{"H", 1, 0, 1},										//Hydrogen's electronegativity is set to 0 because it's not expect to gain or lose electrons
	{"Li", 1, 1.0, 3},										//in a charged diatomic
	{"Be", 2, 1.6, 4},										//have set the default dist to 0.0
	{"Bo", 3, 2.0, 5},										//the first line of the input file should change all of these
	{"C", 4, 2.5, 6},
	{"N", 5, 3.0, 7},
	{"O", 6, 3.5, 8},
	{"F", 7, 4.0, 9},
	{"Na", 1, 0.9, 11},
	{"Mg", 2, 1.3, 12},
	{"Al", 3, 1.6, 13},
	{"Si", 4, 1.9, 14},
	{"P", 5, 2.2, 15},
	{"S", 6, 2.5, 16},
	{"Cl", 7, 3.0, 17}
};

typedef struct {
	char name[5];
	int exp_1;
	int exp_2;
	int exp_3;
	int e_bond;
} mol_comp;

mol_comp comp_ary[82] = {							//table for keeping track of experimental bond order
	{"Al2", 3, 1, 5, 1},
	{"AlCl", 1, 3, 5, 1},
	{"AlF", 1, 3, 5, 1},
	{"AlH", 1, 3, 5, 1},
	{"AlH+", 2, 4, 0, 1},
	{"AlN", 3, 1, 5, 3},
	{"AlO", 2, 4, 6, 1},
	{"AlS", 2, 4, 6, 1},
	{"B2", 3, 1, 5, 1},
	{"BC", 4, 2, 6, 1},
	{"BCl", 1, 3, 5, 1},
	{"BF", 1, 3, 5, 1},
	{"BH", 1, 3, 5, 1},
	{"BN", 3, 1, 5, 3},
	{"BO", 2, 4, 6, 1},
	{"BS", 2, 4, 6, 1},
	{"C2", 1, 3, 5, 2},
	{"CCl", 2, 4, 6, 1},
	{"CF", 2, 4, 6, 1},
	{"CH", 2, 4, 6, 1},
	{"CH-", 3, 1, 5, 1},
	{"Cl2", 1, 3, 5, 1},
	{"Cl2+", 2, 4, 6, 1},
	{"ClF", 1, 3, 5, 1},
	{"ClO", 2, 4, 6, 1},
	{"CN", 2, 4, 6, 3},
	{"CN-", 1, 3, 5, 3},
	{"CO", 1, 3, 5, 3},
	{"CO+", 2, 4, 6, 2},
	{"CP", 2, 4, 6, 3},
	{"CS", 1, 3, 5, 2},
	{"F2", 1, 3, 5, 1},
	{"FO", 2, 4, 6, 1},
	{"H2", 1, 3, 0, 1},
	{"HCl", 1, 3, 5, 1},
	{"HCl+", 2, 4, 6, 1},
	{"HF", 1, 3, 5, 1},
	{"HF+", 2, 4, 6, 1},
	{"N2", 1, 3, 5, 3},
	{"NCl", 3, 1, 5, 2},
	{"NF", 3, 1, 5, 1},
	{"NH", 3, 1, 5, 1},
	{"NH+", 2, 4, 6, 1},
	{"NO", 2, 4, 6, 2},
	{"NO-", 3, 1, 5, 2},
	{"NO+", 1, 3, 5, 3},
	{"NS", 2, 4, 6, 2},
	{"NS+", 1, 3, 5, 3},
	{"O2", 3, 1, 5, 2},
	{"O2+", 2, 4, 6, 2},
	{"OH", 2, 4, 6, 1},
	{"OH-", 1, 3, 5, 1},
	{"OH+", 3, 1, 5, 1},
	{"P2", 1, 3, 5, 3},
	{"P2+", 2, 4, 6, 3},
	{"PCl", 3, 1, 5, 1},
	{"PF", 3, 1, 5, 1},
	{"PF+", 2, 4, 6, 1},
	{"PH", 3, 1, 5, 1},
	{"PH-", 2, 4, 6, 1},
	{"PH+", 2, 4, 6, 1},
	{"PN", 1, 3, 5, 3},
	{"PO", 2, 4, 6, 2},
	{"PO-", 3, 1, 5, 2},
	{"PS", 2, 4, 6, 2},
	{"S2", 3, 1, 5, 2},
	{"S2+", 2, 4, 6, 2},
	{"SCl", 2, 4, 6, 1},
	{"SF", 2, 4, 6, 1},
	{"SH", 2, 4, 6, 1},
	{"SH+", 3, 1, 5, 1},
	{"Si2", 3, 1, 5, 2},
	{"SiCl", 2, 4, 6, 1},
	{"SiF", 2, 4, 6, 1},
	{"SiH", 2, 4, 6, 1},
	{"SiH-", 3, 1, 5, 1},
	{"SiH+", 1, 3, 5, 1},
	{"SiN", 2, 4, 6, 1},
	{"SiO", 1, 3, 5, 1},
	{"SiS", 1, 3, 5, 2},
	{"SO", 3, 1, 5, 2},
	{"SO+", 2, 4, 6, 2}
};