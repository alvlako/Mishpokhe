//#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
//#include <algorithm> // for std::max

float max(float a, float b) {
    return (a > b) ? a : b;
}

int find_max_index(float S_list[], int L_space) {
    int max_index = 0;
    float max_value = S_list[0];
    for (int i = 1; i < L_space; i++) {
        if (S_list[i] > max_value) {
            max_value = S_list[i];
            max_index = i;
        }
    }
    //printf("final_max_index %d\n", max_index);
    return max_index;
}


int find_zero_index_min(float S_list[],  int L_space, int i) {
    int i_1 = 0;
    for (int j = i + 1; j < L_space; j++) {
        if (S_list[j] == 0) {
            i_1 = j;
            break;
        }
    }
    //printf("i_1 in func %d\n", i_1);
    return i_1;
}

int find_zero_index_max(float S_list[],  int L_space, int i, int32_t genome_change_array_c[], int array_space_counter) {
    int start = 0;
    if (array_space_counter == 0){
        int start = -1;
    }
    for (int j = i - 1; j >= 0; j--) {
        if (genome_change_array_c[j] == 1){
            if (S_list[j] == 0) {
                start = j;
            }
            else{
                start = j-1;
            }
            //printf("genome change start j in func %d\n", j);
            break;
        }
        if (S_list[j] == 0) {
            start = j;
            //printf("start j in func %d\n", j);
            break;
        }
    }
    //printf("start in func %d\n", start);
    return start;
}

// Function to compute scores
void ComputeScores(float S_list[], int32_t starts_list[], int32_t i_0, int i_1, float s_x_i_array[], int32_t genome_change_array_c[], int max_gaps) {
    float S_i_min_1 = 0.0f;
    for (int i = i_0; i <= i_1; i++) {
        float s_x_i_flip = s_x_i_array[i];
        if (genome_change_array_c[i] == 1){
            S_i_min_1 = 0;
        }
        float S_i = max(0, S_i_min_1 + s_x_i_flip);
        S_list[i] = S_i;
        S_i_min_1 = S_i;
    }
}


// Function to find cluster matches
int find_cluster_matches(int i_0, int i_1, int k, float S_list[], int starts_list[], int L_space, float S_min, float *cluster_matches_unsorted_s, int *cluster_matches_unsorted_starts, int *cluster_matches_unsorted_ends, float s_x_i_array[], int32_t genome_change_array_c[], int max_gaps) {
    int array_space_counter = 0;
    int i = -1;
    while (1) {
        //printf("i for cumputeScores %d\n", i);
        //printf("i_1 for cumputeScores %d\n",i_1);
        //if (i+1 == i_1){
        //    break;
        //}
        // That is needed to not recompute when it is in another genome
        if (genome_change_array_c[i+1] == 0){
            //printf("non-recomputing ");
            ComputeScores(S_list, starts_list, i+1, i_1, s_x_i_array, genome_change_array_c, max_gaps);
        }
        //printf("S_list array: ");
        //for (int i = 0; i < L_space; i++) {
        //    printf("%f ", S_list[i]);
        //}
        //printf("\n");
        i = find_max_index(S_list, L_space);
        float S_i = S_list[i];
        //printf("S_i %f\n", S_i);
        //printf("S_min %f\n", S_min);
        if (S_i < S_min) {
            break;
        }
        int start = find_zero_index_max(S_list, L_space, i, genome_change_array_c, array_space_counter) +1; //+1
        // that's for the situation when the very first protein has a positive score
        if (start == 1 && S_list[0] > 0){
            //printf("start resolution: ");
            start = start -1;
        }
        cluster_matches_unsorted_s[array_space_counter] = S_i;
        cluster_matches_unsorted_starts[array_space_counter] = start;
        cluster_matches_unsorted_ends[array_space_counter] = i;
        array_space_counter = array_space_counter + 1;
        for (int j = start; j <= i; j++) {
            S_list[j] = 0;
        }
        i_1 = find_zero_index_min(S_list, L_space, i);
        //printf("(%d, %d, %d)\n", cluster_matches_unsorted_s[array_space_counter],  cluster_matches_unsorted_starts[array_space_counter], cluster_matches_unsorted_ends[array_space_counter]);
        //if (array_space_counter > 10){
        //   break;
        //}
    }
    return array_space_counter;
}



