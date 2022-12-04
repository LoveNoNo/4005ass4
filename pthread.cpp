#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <pthread.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"


int size; // problem size
int n_thd; // number of threads
int rowSize; // local number of row of each thread

float* data_odd;
float* data_even;
bool* fire_area;


typedef struct {
    //TODO: specify your arguments for threads
    int thd_id;
    int count;

    //TODO END
} Args;


void initialize(float *data) {
    // intialize the temperature distribution
    int len = size * size;
    for (int i = 0; i < len; i++) {
        data[i] = wall_temp;
    }
}


void generate_fire_area(bool *fire_area){
    // generate the fire area
    int len = size * size;
    for (int i = 0; i < len; i++) {
        fire_area[i] = 0;
    }

    float fire1_r2 = fire_size * fire_size;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            int a = i - size / 2;
            int b = j - size / 2;
            int r2 = 0.5 * a * a + 0.8 * b * b - 0.5 * a * b;
            if (r2 < fire1_r2) fire_area[i * size + j] = 1;
        }
    }

    float fire2_r2 = (fire_size / 2) * (fire_size / 2);
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            int a = i - 1 * size / 3;
            int b = j - 1 * size / 3;
            int r2 = a * a + b * b;
            if (r2 < fire2_r2) fire_area[i * size + j] = 1;
        }
    }
}


void update(float *data, float *new_data, int tid) {
    // TODO: update the temperature of each point, and store the result in `new_data` to avoid data racing
    int startRow = tid*(size-2)/n_thd + 1;

    for (int i = startRow; i < rowSize -1; i++){ // 计算范围：up most & down most not write result
        for (int j = 1; j < size - 1; j++){  // 计算顺序：一维数组，逐行扫描
            int idx = i * size + j;
            float up = data[idx - size];
            float down = data[idx + size];
            float left = data[idx - 1];
            float right = data[idx + 1];
            float new_val = (up + down + left + right) / 4; // 计算方式：当前点温度等于上下左右四个点温度的平均值
            new_data[idx] = new_val;
        }
    }
}


void maintain_fire(float *data, bool* fire_area) {
    // TODO: maintain the temperature of fire
    int len = rowSize * size;
    for (int i = 0; i < len; i++){
        if (fire_area[i]) data[i] = fire_temp; 
    }
}


void maintain_wall(float *data, int tid) {
    // TODO: maintain the temperature of the wall
    int startRow = tid*(size-2)/n_thd + 1;

    for(int i=0; i<size; i++){ // only first/last thread maintain first/last row
        if(tid == 0){
            data[i] = wall_temp; // 第一行 
        }
        else if(tid == n_thd-1){ // only last thread maintain last row 
            data[(size-1)*size + i] = wall_temp; //最后一行
        }   
    }

    for(int j=startRow; j<rowSize-1; j++){
        data[j*size] = wall_temp;   // 最左列
        data[j*size + (size -1)] = wall_temp; //最右列
    }
}

void* worker(void* args) {
    // Pthread TODO: procedure in each threads
    Args* my_arg = (Args*) args;
    int tid = my_arg->thd_id;
    int count = my_arg->count;

    if(count % 2 ==1){
        update(data_odd, data_even, tid);
        maintain_fire(data_even, fire_area);
        maintain_wall(data_even, tid);
    }else{
        update(data_even, data_odd, tid);
        maintain_fire(data_odd, fire_area);
        maintain_wall(data_odd, tid);
    }

    // TODO END
}

int main(int argc, char *argv[]) {
    size = atoi(argv[1]);
    n_thd = atoi(argv[2]);
    int rowSize = (size-2)/n_thd +2;

    pthread_t thds[n_thd];

    Args args[n_thd];
    for (int thd = 0; thd < n_thd; thd++){
        args[thd].thd_id = thd;
    }

    data_odd = new float[size * size];
    data_even = new float[size * size];
    fire_area = new bool[size * size];

    generate_fire_area(fire_area);
    initialize(data_odd);

    int count = 1;
    double total_time = 0;
    

    while (true) {
        for (int thd = 0; thd < n_thd; thd++){
            args[thd].count = count;
        }

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        
        // PThread TODO: assign jobs
        
        for (int thd = 0; thd < n_thd; thd++) pthread_create(&thds[thd], NULL, worker, &args[thd]);
        for (int thd = 0; thd < n_thd; thd++) pthread_join(thds[thd], NULL);

        //TODO End

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        double this_time = std::chrono::duration<double>(t2 - t1).count();
        total_time += this_time;
        printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        count++;

        if(count > 1000){break;}

    }

    printf("Converge after %d iterations, elapsed time: %.6f, average computation time: %.6f\n", count-1, total_time, (double) total_time / (count-1));
    
    printf("Student ID: 117010332\n"); // replace it with your student id
    printf("Name: XU Jiale\n"); // replace it with your name
    printf("Assignment 4: Heat Distribution Simulation Pthread Implementation\n");

    delete[] data_odd;
    delete[] data_even;
    delete[] fire_area;

	return 0;
}