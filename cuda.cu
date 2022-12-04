#include <cuda.h>
#include <cuda_runtime.h>

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


int block_size = 512; // cuda thread block size
int size; // problem size
int n_thd = 4; // number of thread


__global__ void initialize(float *data) {
    // TODO: intialize the temperature distribution (in parallelized way)
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    data[i] = wall_temp;
}


__global__ void generate_fire_area(bool *fire_area){
    // TODO: generate the fire area (in parallelized way)
    
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


__global__ void update(float *data, float *new_data, int size) {
    // TODO: update temperature for each point  (in parallelized way)
    
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    
    if(i > size && i < size*size-size){ //第一行和最后一行的上下会超过data范围
        float up = data[i - size];
        float down = data[i + size];
        float left = data[i - 1];
        float right = data[i + 1];
        float new_val = (up + down + left + right) / 4; 
        new_data[i] = new_val;
    }
}


__global__ void maintain_wall(float *data, int size) {
    // TODO: maintain the temperature of the wall (sequential is enough)

    for(int i=0; i<size; i++){
        data[i] = wall_temp; // 第一行 
        data[(size-1)*size + i] = wall_temp; //最后一行
        data[i*size] = wall_temp;   // 最左列
        data[i*size + (size -1)] = wall_temp; //最右列
    }
}


__global__ void maintain_fire(float *data, bool *fire_area) {
    // TODO: maintain the temperature of the fire (in parallelized way)
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    
    if(fire_area[i]) data[i] = fire_temp;
}


#ifdef GUI
__global__ void data2pixels(float *data, GLubyte* pixels){
    // TODO: convert rawdata (large, size^2) to pixels (small, resolution^2) for faster rendering speed (in parallelized way)
}


void plot(GLubyte* pixels){
    // visualize temprature distribution
    #ifdef GUI
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(resolution, resolution, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    glutSwapBuffers();
    #endif
}
#endif


void master() {
    float *data_odd;
    float *data_even;
    bool *fire_area;

    cudaMalloc(&data_odd, size * size * sizeof(float));
    cudaMalloc(&data_even, size * size * sizeof(float));
    cudaMalloc(&fire_area, size * size * sizeof(bool));

    #ifdef GUI
    GLubyte *pixels;
    GLubyte *host_pixels;
    host_pixels = new GLubyte[resolution * resolution * 3];
    cudaMalloc(&pixels, resolution * resolution * 3 * sizeof(GLubyte));
    #endif

    int n_block_size = size * size / block_size + 1;
    //int n_block_resolution = resolution * resolution / block_size + 1;

    initialize<<<n_block_size, block_size>>>(data_odd);
    generate_fire_area<<<1, 1>>>(fire_area);

    // initialize<<<1, n_thd>>>(data_odd);
    // generate_fire_area<<<1, n_thd>>>(fire_area);
    
    int count = 1;
    double total_time = 0;

    while (true){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // TODO: modify the following lines to fit your need.
        
        if (count % 2 == 1) {
            update<<<n_block_size, block_size>>>(data_odd, data_even, size);
            maintain_fire<<<n_block_size, block_size>>>(data_even, fire_area);
            maintain_wall<<<1, 1>>>(data_even, size);
        } else {
            update<<<n_block_size, block_size>>>(data_even, data_odd, size);
            maintain_fire<<<n_block_size, block_size>>>(data_odd, fire_area);
            maintain_wall<<<1, 1>>>(data_odd, size);
        }

        /* One block setting thread
        if (count % 2 == 1) {
            update<<<1, n_thd>>>(data_odd, data_even);
            maintain_fire<<<1, n_thd>>>(data_even, fire_area);
            maintain_wall<<<1, 1>>>(data_even);
        } else {
            update<<<1, n_thd>>>(data_even, data_odd);
            maintain_fire<<<1, n_thd>>>(data_odd, fire_area);
            maintain_wall<<<1, 1>>>(data_odd);
        }*/

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        double this_time = std::chrono::duration<double>(t2 - t1).count();
        total_time += this_time;
        printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        count++;
        
        #ifdef GUI
        if (count % 2 == 1) {
            data2pixels<<<n_block_resolution, block_size>>>(data_even, pixels);
        } else {
            data2pixels<<<n_block_resolution, block_size>>>(data_odd, pixels);
        }
        cudaMemcpy(host_pixels, pixels, resolution * resolution * 3 * sizeof(GLubyte), cudaMemcpyDeviceToHost);
        plot(host_pixels);
        #endif

    }

    printf("Converge after %d iterations, elapsed time: %.6f, average computation time: %.6f\n", count-1, total_time, (double) total_time / (count-1));


    cudaFree(data_odd);
    cudaFree(data_even);
    cudaFree(fire_area);

    #ifdef GUI
    cudaFree(pixels);
    delete[] host_pixels;
    #endif
    
}


int main(int argc, char *argv[]){
    
    size = atoi(argv[1]);
    n_thd = atoi(argv[2]);

    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(resolution, resolution);
    glutCreateWindow("Heat Distribution Simulation Sequential Implementation");
    gluOrtho2D(0, resolution, 0, resolution);
    #endif

    master();

    printf("Student ID: 117010332\n"); // replace it with your student id
    printf("Name: XU Jiale\n"); // replace it with your name
    printf("Assignment 4: Heat Distribution CUDA Implementation\n");

    return 0;

}


