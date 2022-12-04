#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"


int size;
int n_omp_threads;

float* data_odd;
float* data_even;
bool* fire_area;


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


void update(float *data, float *new_data) {
    // TODO: update the temperature of each point, and store the result in `new_data` to avoid data racing
    
    #pragma omp parallel for private(j) firstprivate(data) lastprivate(new_data)
    for (int i = 1; i < size - 1; i++){ 
        for (int j = 1; j < size - 1; j++){  
            int idx = i * size + j;
            float up = data[idx - size];
            float down = data[idx + size];
            float left = data[idx - 1];
            float right = data[idx + 1];
            float new_val = (up + down + left + right) / 4;
            new_data[idx] = new_val;
        }
    }
}


void maintain_fire(float *data, bool* fire_area) {
    // TODO: maintain the temperature of fire
    int len = size * size;

    #pragma omp parallel for firstprivate(fire_area) lastprivate(data)
    for (int i = 0; i < len; i++){
        if (fire_area[i]) data[i] = fire_temp; 
    }

}


void maintain_wall(float *data) {
    // TODO: maintain the temperature of the wall
    for(int i=0; i<size; i++){
        data[i] = wall_temp; 
        data[(size-1)*size + i] = wall_temp; 
        data[i*size] = wall_temp;   
        data[i*size + (size -1)] = wall_temp; 
    }

}


int main(int argc, char *argv[]){
    
    size = atoi(argv[1]);
    n_omp_threads = atoi(argv[2]);

    data_odd = new float[size * size];
    data_even = new float[size * size];
    fire_area = new bool[size * size];

    generate_fire_area(fire_area);
    initialize(data_odd);

    int count = 1;
    double total_time = 0;   

    while (true) {

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        
        if (count % 2 == 1) {
            update(data_odd, data_even);
            maintain_fire(data_even, fire_area); // 炉火和墙体的温度在update计算后强制刷新为固定温度
            maintain_wall(data_even);
        } else {
            update(data_even, data_odd);
            maintain_fire(data_odd, fire_area);
            maintain_wall(data_odd);
        }
        
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
    printf("Assignment 4: Heat Distribution Simulation OpenMP Implementation\n");

    delete[] data_odd;
    delete[] data_even;
    delete[] fire_area;

	return 0;
}
