#include "exp_data.h"
#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <math.h>

#define MAX 5000
char folder_name[100];
int list_size = 0;
int isStable = 0;

typedef struct{
	double x;
	double y;
	double z;
} coor;

typedef struct{
	char name[100];
	char file_name[200];
	char mol_name[6];
	double energy;
	int bond;
	int mult;
} structure;
structure mins[MAX];
structure srt_mins[MAX];

double get_coor(char string[7]){
	char substr[6];
	if(string[0] != '-'){
		memcpy(substr, &string[1], 5);
		substr[5] = '\0';
		return strtod(substr, NULL);
	}
	return strtod(string, NULL);
}

float calc_bond(char f_name[200]){
	coor mol[50];
	float bond = 0;
	int cur = 0;
	double len;
	char buf[100];
	char subbuf[7];
	char out[10];
	int start = 31;
	FILE *f_ground = fopen(f_name, "r");
	fgets(buf, 100, f_ground);
	while (buf[0] != 'H') {
		fgets(buf, 100, f_ground);
	}
	fgets(buf, 100, f_ground);
	fgets(buf, 100, f_ground);
	while(strcmp(buf, "END\n") != 0){
		memcpy(subbuf, &buf[31], 6);
		subbuf[6] = '\0';
		mol[cur].x = get_coor(subbuf);
	
		memcpy(subbuf, &buf[39], 6);
		subbuf[6] = '\0';
		mol[cur].y = get_coor(subbuf);
	
		memcpy(subbuf, &buf[47], 6);	
		subbuf[6] = '\0';
		mol[cur].z = get_coor(subbuf);
		cur++;
		fgets(buf, 100, f_ground);
	}
	
	len = fabs(mol[0].z);
	if (len < fabs(mol[1].z)) {
		len = fabs(mol[1].z);
	}
	len *= 1.5;
	for(int i = 2; i < cur; i++){
		if(fabs(mol[i].x) < len && fabs(mol[i].y) < len){
			 if ((0.8 * mol[1].z < mol[i].z && mol[i].z < 0.8 * mol[0].z) || (0.8 * mol[0].z < mol[i].z && mol[i].z < 0.8 * mol[1].z)) {
			 	bond++;
			 }
		}
		
		if(fabs(mol[i].x) < len && fabs(mol[i].y) < len && isStable == 0){
			 if (!((2 * mol[1].z < mol[i].z && mol[i].z < 2 * mol[0].z) || (2 * mol[0].z < mol[i].z && mol[i].z < 2 * mol[1].z))) {
			 	isStable = -1;
			 } 
		} else {
			isStable = -1;
		}
	}
	fclose(f_ground);
	bond /= 2.0;
	return bond;
}
/*
void sort_mins(){
	for(int i = 0; i < list_size; i++) {
		if(mins[i].){
		
		}
	}
}*/

void output(){
	//sort_mins();
	char buf[1000];
	char tbl_name[200];
	int unstable = 0;
	int ord_er = 0;
	int bond_er = 0;
	char out_str[169];
	char empt_str[169];
	int ld_ord[3] = {-1,-1,-1};
	int ld_mults[3] = {0,0,0};
	int num_mults = 0;
	int j = 0;
	int k = 0;
	int ary_ind = 0;
	char temp_str[50];
	char new_name[200];
	char f_name[200];
	char g_name[100];
	char bond[10];
	FILE *Source;
	FILE *Dest;
	
	memset(empt_str, ' ', 169);
	empt_str[9] = '|';
	empt_str[10] = '|';
	empt_str[23] = '|';
	empt_str[24] = '|';
	empt_str[42] = '|';
	empt_str[43] = '|';
	empt_str[55] = '|';
	empt_str[56] = '|';
	empt_str[73] = '|';
	empt_str[74] = '|';
	empt_str[88] = '|';
	empt_str[89] = '|';
	empt_str[102] = '|';
	empt_str[103] = '|';
	empt_str[114] = '|';
	empt_str[115] = '|';
	empt_str[132] = '|';
	empt_str[133] = '|';
	empt_str[150] = '|';
	empt_str[151] = '|';
	empt_str[167] = '\n';
	empt_str[168] = '\0';
	
	strcpy(g_name, "Min_Structs_");
	strcat(g_name, folder_name);	
	strcat(g_name, "/Grounds");
	if(access(g_name, F_OK) == -1) {
		mkdir(g_name, 0777);
	} else {
		DIR *MinDir;
		MinDir = opendir(g_name);
	}
	
	strcpy(tbl_name, "Min_Structs_");
	strcat(tbl_name, folder_name);
	strcat(tbl_name, "/Results_Table_");
	strcat(tbl_name, folder_name);
	strcat(tbl_name, ".txt");
	FILE *output = fopen(tbl_name, "w");
	fprintf(output, "Mol Name || Exp. Order || Lewis Dot Order || Exp. Bond || Lewis Dot Bond || Order Error || Bond Error ||  Stable  ||   Min Energy   ||   Mid Energy   ||   Max Energy  \n");
	fprintf(output, "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	for(int i = 0; i < list_size; i++){
		while(strcmp(mins[i].mol_name, comp_ary[ary_ind].name) != 0) {
			ary_ind++;
		}
		strncpy(out_str, empt_str, 169);
		
		while(mins[i].mol_name[j] != '\0') {
			out_str[j+2] = mins[i].mol_name[j];
			j++;
		}
		
		j = 14;
		if(comp_ary[ary_ind].exp_3 != 0) {
			out_str[j++] = comp_ary[ary_ind].exp_1 + '0';
			out_str[j++] = ',';
			out_str[j++] = comp_ary[ary_ind].exp_2 + '0';
			out_str[j++] = ',';
			out_str[j++] = comp_ary[ary_ind].exp_3 + '0';
		} else {
			out_str[j++] = comp_ary[ary_ind].exp_1 + '0';
			out_str[j++] = ',';
			out_str[j++] = comp_ary[ary_ind].exp_2 + '0';
		}
		
		ld_ord[0] = i;
		if(strcmp(mins[i].mol_name, mins[i+1].mol_name) == 0){
			num_mults++;
			if(mins[i].energy < mins[i+1].energy) {
				ld_ord[1] = i+1;
			} else {
				ld_ord[0] = i+1;
				ld_ord[1] = i;
			}
		}
		
		if(strcmp(mins[i].mol_name, mins[i+2].mol_name) == 0){
			num_mults++;
			if(mins[ld_ord[1]].energy < mins[i+2].energy) {
				ld_ord[2] = i+2;
			} else if(mins[ld_ord[0]].energy < mins[i+2].energy) {
				ld_ord[2] = ld_ord[1];
				ld_ord[1] = i+2;
			} else {
				ld_ord[2] = ld_ord[1];
				ld_ord[1] = ld_ord[0];
				ld_ord[0] = i+2;
			}
		}
		
		ld_mults[0] = mins[ld_ord[0]].mult;
		ld_mults[1] = mins[ld_ord[1]].mult;
		ld_mults[2] = mins[ld_ord[2]].mult;
		
		j = 31;
		out_str[j++] = mins[ld_ord[0]].mult + '0';
		if(mins[ld_ord[1]].mult != 0) {
			out_str[j++] = ',';
			out_str[j++] = mins[ld_ord[1]].mult + '0';
		}
		if(mins[ld_ord[2]].mult != 0 && comp_ary[ary_ind].exp_3 != 0) {
			out_str[j++] = ',';
			out_str[j++] = mins[ld_ord[2]].mult + '0';
		}
		
		memset(new_name, '\0', 200);
		memset(f_name, '\0', 200);
		strncpy(new_name, mins[ld_ord[0]].file_name, 200);
		new_name[strlen(new_name)-11] = '\0';
		strncat(new_name, "ground_MCfinal.pdb", 200 - strlen(new_name));
		rename(mins[ld_ord[0]].file_name, new_name);
		
		Source = fopen(new_name, "r");
		
		k = 0;
		strcpy(f_name, "Min_Structs_");
		strcat(f_name, folder_name);
		strcat(f_name, "/Grounds/");	
		while(new_name[k] != '/'){
			k++;
		}
		k++;
		while(new_name[k] != '\0'){
			f_name[strlen(f_name)] = new_name[k];
			k++;
		}
		Dest = fopen(f_name, "w+");
		while (fgets(buf, sizeof(buf), Source)) {
			fprintf(Dest, "%s", buf);
		}
		fclose(Source);
		fclose(Dest);
		
		out_str[49] = comp_ary[ary_ind].e_bond + '0';
		
		sprintf(bond, "%.1f", calc_bond(new_name));
		out_str[64] = bond[0];
		if(bond[2] != '0'){
			out_str[65] = bond[1];
			out_str[66] = bond[2];
		}
		
		if(out_str[49] != out_str[64] || out_str[50] != out_str[65]){
			out_str[95] = '1';
			bond_er++;
		} else {
			out_str[95] = '0';
		}
		
		j = 81;
		if ((comp_ary[ary_ind].exp_1 - ld_mults[0]) != 0 ||
			(comp_ary[ary_ind].exp_2 - ld_mults[1]) != 0 ||
			(((comp_ary[ary_ind].exp_3 - ld_mults[2]) != 0) && 
			(comp_ary[ary_ind].exp_3 != 0))){
			out_str[j] = '1';
			ord_er++;
		} else {
			out_str[j] = '0';
		}
		
		if(isStable == -1){
			out_str[108] = '1';
			unstable++;
		} else {
			out_str[108] = '0';
		}
		
		j = 117;
		k = 0;
		snprintf(temp_str, 20, "%f", mins[ld_ord[0]].energy);
		while(temp_str[k] != '\0') {
			out_str[j] = temp_str[k];
			j++;
			k++;
		}
		
		if (ld_ord[1] != -1) {
			j = 135;
			k = 0;
			snprintf(temp_str, 20, "%f", mins[ld_ord[1]].energy);
			while(temp_str[k] != '\0') {
				out_str[j] = temp_str[k];
				j++;
				k++;
			}
		}
		
		if (ld_ord[2] != -1 && comp_ary[ary_ind].exp_3 != 0) {
			j = 153;
			k = 0;
			snprintf(temp_str, 20, "%f", mins[ld_ord[2]].energy);
			while(temp_str[k] != '\0') {
				out_str[j] = temp_str[k];
				j++;
				k++;
			}
		}
		
		fprintf(output, "%s", out_str);
		j = 0;
		k = 0;
		memset(ld_ord, -1, 3);
		memset(ld_mults, 0, 3);
		i += num_mults;
		num_mults = 0;
		isStable = 0;
	}
	fprintf(output, "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	strncpy(out_str, empt_str, 169);
	j = 81;
	k = 0;
	strncpy(out_str, "Total    ||     -      ||        -        ||     -     ||        -       ||             ||            ||          ||       -        ||       -        ||       -       \n", 169);
	snprintf(temp_str, 4, "%d", ord_er);
	while(temp_str[k] != '\0') {
		out_str[j] = temp_str[k];
		j++;
		k++;
	}
	j = 95;
	k = 0;
	snprintf(temp_str, 4, "%d", bond_er);
	while(temp_str[k] != '\0') {
		out_str[j] = temp_str[k];
		j++;
		k++;
	}
	j = 108;
	k = 0;
	snprintf(temp_str, 4, "%d", unstable);
	while(temp_str[k] != '\0') {
		out_str[j] = temp_str[k];
		j++;
		k++;
	}
	fprintf(output, "%s", out_str);
	fclose(output);
}

void copy_mins(){
	int under_count = 0;
	int j = 0;
	char buf[1000];
	char f_name[200];
	char new_name[200];
	char fold_name[200];
	FILE *Source;
	FILE *Dest;
	
	strcpy(fold_name, "Min_Structs_");
	strcat(fold_name, folder_name);
	if(access(fold_name, F_OK) == -1) {
		mkdir(fold_name, 0777);
	} else {
		DIR *MinDir;
		MinDir = opendir(fold_name);
	}
	for(int i = 0; i < list_size; i++) {
		while(under_count < 1){
			new_name[j] = mins[i].name[j];
			if (mins[i].name[j] == '_') {/////////////////////////////////////////////
				under_count++;
			}
			j++;
		}
		if (mins[i].name[j] == 'p' || mins[i].name[j] == 'm') {
			new_name[j] = mins[i].name[j];
			new_name[j+1] = mins[i].name[j+1];
			j += 2;
		}
		
		new_name[j] = mins[i].name[j];
		new_name[j+1] = mins[i].name[j+1];
		new_name[j+2] = mins[i].name[j+2];
		new_name[j+3] = mins[i].name[j+3];
		new_name[j+4] = mins[i].name[j+4];
		
		mins[i].mult += (new_name[j+0] - '0') - (new_name[j+2] - '0') + 1;
		strcat(new_name, "min_MCfinal.pdb");
		j = 0;
		under_count = 0;
		
		strncpy(f_name, folder_name, 200);
		strcat(f_name, "/");
		strcat(f_name, mins[i].name);
		Source = fopen(f_name, "r");
	
		strncpy(f_name, fold_name, 200);
		strncat(f_name, "/", 100);
		strcat(f_name, new_name);
		Dest = fopen(f_name, "w+");
		strcpy(mins[i].file_name, f_name);
		while (fgets(buf, sizeof(buf), Source)) {
			fprintf(Dest, "%s", buf);
		}
		
		memset(new_name, '\0', 200);
		fclose(Source);
		fclose(Dest);
	}	
}

void comp_energies(char Struct_array[25][100]){
	char buf[100];
	int i = 0;
	int j = 25;
	int k = 4;
	char min_struct[100];
	double min_energy = 1;
	double cur_energy;
	char energy[100];
	char cur_struct[200];
	FILE *file;
	
	for (int i = 0; i < 25; i++) {
		j = 25;
		memset(buf, '\0', 100);
		if(Struct_array[i][0] == 's'){
			strncpy(cur_struct, folder_name, sizeof(folder_name));
			strcat(cur_struct, "/");
			strncat(cur_struct, Struct_array[i], sizeof(Struct_array[i]));
			file = fopen(cur_struct, "r");
			while (buf[0] != 'H') {
				fgets(buf, 100, file);
			}
			int ln = strlen(buf) - 1;
			if (buf[ln] == '\n') {
				buf[ln] = '\0';
			}
			while (buf[j] != ' ' && buf[j] != '\n'&& buf[j] != '\0'){
				energy[j-25] = buf[j];
				j++;
			}
			cur_energy = strtod(energy, NULL);
			if (min_energy > cur_energy || min_energy == 1){
				min_energy = cur_energy;
				strncpy(min_struct, Struct_array[i], sizeof(Struct_array[i]));
			}
		}
		fclose(file); 	
	}
	
	while(min_struct[k] != '_') {///////////////////////////////////////////////////////////////
        mins[list_size].mol_name[k-4] = min_struct[k];
        k++;
    }
    
    if(min_struct[k+1] == 'p'){
        mins[list_size].mol_name[k-4] = '+';
	} else if(min_struct[k+1] == 'm'){
        mins[list_size].mol_name[k-4] = '-';
	}
	
    strcpy(mins[list_size].name, min_struct);
	mins[list_size].energy = min_energy;
	list_size++;
}

void analysis(char Struct_names[MAX][100]) {
	int count = 0;
	int current = 0;
	int under_count = 0;
	int i = 0;
	int j = 0;
	int mult;
	int done_once = 0;
	int charged = 0;
	char cur_string[100];
	char short_str[100];
	char Struct_array[25][100];
	while (count < MAX && strcmp(Struct_names[count], "") != 0){
		if(	Struct_names[count][0] == 's' &&
			Struct_names[count][1] == 'p' &&
			Struct_names[count][2] == 'i' &&
			Struct_names[count][3] == 'n'){
			while(under_count < 1){
				short_str[i] = Struct_names[count][j];
				if (Struct_names[count][j] == '_') {////////////////////////////////////////
					under_count++;
				}
				i++;
				j++;
			}
			if (Struct_names[count][j] == 'p' || Struct_names[count][j] == 'm') {
				short_str[i] = Struct_names[count][j];
				short_str[i+1] = Struct_names[count][j+1];
				i += 2;
				j += 2;
			}
			
			short_str[i] = Struct_names[count][j];
			short_str[i+1] = Struct_names[count][j+1];
			short_str[i+2] = Struct_names[count][j+2];
			short_str[i+3] = Struct_names[count][j+3];
			short_str[i+4] = Struct_names[count][j+4];
			
			i = 0;
			j = 0;
			under_count = 0;
			if(strcmp(cur_string, short_str) == 0){
				strcpy(Struct_array[current], Struct_names[count]);
				current++;
			} else {
				if (done_once != 0) {
					comp_energies(Struct_array);
				} else {
					done_once++;
				}
				for(int k = 0; k < 25; k++){
					memset(Struct_array[k], '\0', 100);
				}
				current = 1;
				strcpy(Struct_array[0], Struct_names[count]);
				strcpy(cur_string, short_str); 
			}
		}
		count++;
	}
	comp_energies(Struct_array);
	copy_mins();
	output();
}

void read_sort(char folder_name[100]){
	struct dirent Struct_dirent[MAX + 2];
	char Struct_names[MAX][100];
	int count = 0;
	int cur = 0;
	DIR *StructsDir;
	StructsDir = opendir(folder_name);
	struct dirent *file_name;
	if (StructsDir != NULL) {
		file_name = readdir(StructsDir);
		while (file_name && (count < MAX + 2)) {
			Struct_dirent[count] = *file_name;
			count++;
			file_name = readdir(StructsDir);
        }
        file_name = readdir(StructsDir);
        if (file_name){
        	printf("Increase MAX, too many files");
        	exit(-1);
        }
		closedir (StructsDir);
	} else {
    	perror ("Couldn't open the directory");
    	exit(-1);
	}
	
	count = 0;
	while (count < MAX){
		strcpy(Struct_names[count], Struct_dirent[count + 2].d_name);
		count++;
	}
	analysis(Struct_names);
}

int main() {
	char buf[100];										//buffer for reading input file
	char quit[] = "q";
	int cmp;
	printf("Final Structures Folder Name (q to quit): \n");		//prompts for the name of the folder
	fgets(buf, 100, stdin);								//takes in user input
	strcpy(folder_name, buf);							//copy the name of the file from the buffer
	int ln = strlen(folder_name) - 1;					//finds end of string and removes the new line
	if (folder_name[ln] == '\n') {
    	folder_name[ln] = '\0';
    }
    cmp = strcmp(folder_name, quit);
	if(cmp == 0){
		printf("\n");
		return 0;
	}
	if(access(folder_name, F_OK) == -1) {
		printf("Folder Not Found\n\n");
		return -1;
	}
	read_sort(folder_name);
	return 0;
}