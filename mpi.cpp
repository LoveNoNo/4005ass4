#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"


int size; // problem size


int my_rank;
int world_size;
int rowSize; // local number of row of each thread


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


void update(float *data, float *new_data/*, int begin, int end*/) {
    // TODO: update the temperature of each point, and store the result in `new_data` to avoid data racing
    for (int i = 1; i < rowSize - 1; i++){ // 计算范围：最外圈一层是wall或者fire 温度一直保持不变，所以从1开始到size-2结束
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


void maintain_wall(float *data) {
    // TODO: maintain the temperature of the wall
    for(int i=0; i<rowSize; i++){ // only first thread maintain first row
        if(my_rank == 0){
            data[i] = wall_temp; // 第一行 
        }
        if(my_rank == world_size-1){ // only last thread maintain last row 
            data[(size-1)*size + i] = wall_temp; //最后一行
        }
        data[i*size] = wall_temp;   // 最左列
        data[i*size + (size -1)] = wall_temp; //最右列
    }
}


#ifdef GUI
void data2pixels(float *data, GLubyte* pixels, int begin, int end){
    // convert rawdata (large, size^2) to pixels (small, resolution^2) for faster rendering speed
    float factor_data_pixel = (float) size / resolution;
    float factor_temp_color = (float) 255 / fire_temp;
    for (int x = 0; x < resolution; x++){
        for (int y = 0; y < resolution; y++){
            int idx = x * resolution + y;
            int idx_pixel = idx * 3;
            int x_raw = x * factor_data_pixel;
            int y_raw = y * factor_data_pixel;
            int idx_raw = y_raw * size + x_raw;
            float temp = data[idx_raw];
            int color =  ((int) temp / 5 * 5) * factor_temp_color;
            pixels[idx_pixel] = color;
            pixels[idx_pixel + 1] = 255 - color;
            pixels[idx_pixel + 2] = 255 - color;
        }
    }
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



void slave(){
    // TODO: MPI routine (one possible solution, you can use another partition method)
    //int my_begin_row_id = size * my_rank / (world_size);
    //int my_end_row_id = size * (my_rank + 1) / world_size;

    
    float* data_odd = new float[size * size];
    bool* fire_area = new bool[size * size];
    int startID = my_rank*(size-2)/world_size*size;

    // TODO: Initialize a storage for temperature distribution of this area
    float* local_data_odd = new float[rowSize * size];
    float* local_data_even = new float[rowSize * size];
    bool* fire_area_local = new bool[rowSize * size];
    float sendBuf_up[size];
    float sendBuf_down[size];
    float recvBuf_up[size];
    float recvBuf_down[size];

    // TODO: Receive initial temperature distribution of this area from master
    MPI_Bcast(data_odd, size*size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(fire_area, size*size, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    // get each local data from whole data array
    for(int i=0; i<rowSize*size; i++){
        local_data_odd[i] = data_odd[startID + i];
        fire_area_local[i] = fire_area[startID + i];
    }

    // TODO: Initialize a storage for local pixels (pls refer to sequential version for initialization of GLubyte)
    #ifdef GUI
    GLubyte* local_pixcels;
    #endif

    bool cont = true;
    int count_local = 1;
    while (cont) {
        // TODO: computation part
        if (count_local % 2 == 1) {
            update(local_data_odd, local_data_even);
            maintain_fire(local_data_odd, fire_area_local);
            maintain_wall(local_data_odd);
        }else{
            update(local_data_even, local_data_odd);
            maintain_fire(local_data_even, fire_area_local);
            maintain_wall(local_data_even);
        }
        count_local++;

        // TODO: after computation, send border row to neighbours
        for(int i=0; i<size; i++){
            if(count_local % 2 ==1){
                sendBuf_up[i] = local_data_odd[i];
                sendBuf_down[i] = local_data_odd[(rowSize -1)*size + i];
            }else{
                sendBuf_up[i] = local_data_even[i];
                sendBuf_down[i] = local_data_even[(rowSize -1)*size + i];
            }
        }

        // recv up & send down, slave thread frist recv
        // recv down & send up; last thread first send
        MPI_Recv(recvBuf_up, size, MPI_FLOAT, my_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if(my_rank < world_size-1){
            MPI_Send(sendBuf_down, size, MPI_FLOAT, my_rank+1, 1, MPI_COMM_WORLD);
            MPI_Recv(recvBuf_down, size, MPI_FLOAT, my_rank+1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(sendBuf_up, size, MPI_FLOAT, my_rank-1, 2, MPI_COMM_WORLD);
        }else{
            MPI_Send(sendBuf_up, size, MPI_FLOAT, my_rank-1, 2, MPI_COMM_WORLD);
        }
        
        for(int i=0; i<size; i++){
            if(count_local % 2 ==1){
                local_data_even[i] = recvBuf_up[i];
                if(my_rank < world_size-1){
                    local_data_even[(rowSize -1)*size + i] = recvBuf_down[i];
                }
            }else{
                local_data_odd[i] = recvBuf_up[i];
                if(my_rank < world_size-1){
                    local_data_odd[(rowSize -1)*size + i] = recvBuf_down[i];
                }
            }
        }
        #ifdef GUI
        // TODO: conver raw temperature to pixels (much smaller than raw data)

        // TODO: send pixels to master (you can use MPI_Byte to transfer anything to master, then you won't need to declare MPI Type :-) )

        #endif

    }

    #ifdef GUI
    data2pixels(local_data, local_pixcels);
    #endif

    // TODO: Remember to delete[] local_data and local_pixcels.
}



void master() {
    // TODO: MPI routine (one possible solution, you can use another partition method)
    float* data_odd = new float[size * size];
    //float* data_even = new float[size * size];
    bool* fire_area = new bool[size * size];
    bool* fire_area_local = new bool[rowSize * size];

    int rowSize = (size-2)/world_size +2;
    int startID = my_rank*(size-2)/world_size*size;

    float* local_data_odd = new float[rowSize * size];
    float* local_data_even = new float[rowSize * size];

    float sendBuf_down[size];
    float recvBuf_down[size];

    initialize(data_odd);
    generate_fire_area(fire_area);

    #ifdef GUI
    GLubyte* pixels;
    pixels = new GLubyte[resolution * resolution * 3];
    #endif

    int count = 1;
    double total_time = 0;

    // TODO: Send initial distribution to each slave process
    MPI_Bcast(data_odd, size*size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(fire_area, size*size, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

    for(int i=0; i<rowSize*size; i++){
        local_data_odd[i] = data_odd[startID + i];
        fire_area_local[i] = fire_area[startID + i];
    }
    
    while (true) {
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // TODO: Computation of my part
        if (count % 2 == 1) {
            update(local_data_odd, local_data_even);
            maintain_fire(local_data_odd, fire_area_local);
            maintain_wall(local_data_odd);
        }else{
            update(local_data_even, local_data_odd);
            maintain_fire(local_data_even, fire_area_local);
            maintain_wall(local_data_even);
        }
        
        // TODO: Send border row to neighbours
        for(int i=0; i<size; i++){
            if(count % 2 ==1){
                sendBuf_down[i] = local_data_odd[(rowSize -1)*size + i];
            }else{
                sendBuf_down[i] = local_data_even[(rowSize -1)*size + i];
            }
        }

        MPI_Send(sendBuf_down, size, MPI_FLOAT, 1, 1, MPI_COMM_WORLD);
        MPI_Recv(recvBuf_down, size, MPI_FLOAT, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for(int i=0; i<size; i++){
            if(count % 2 ==1){
                local_data_even[(rowSize -1)*size + i] = recvBuf_down[i];
            }else{
                local_data_odd[(rowSize -1)*size + i] = recvBuf_down[i];
            }
        }

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        double this_time = std::chrono::duration<double>(t2 - t1).count();
        total_time += this_time;
        printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        count++;

        if(count > 1000){break;}

        #ifdef GUI
        if (count % 2 == 1) {
            // TODO: Gather pixels of slave processes
            data2pixels(data_even, pixels);
        } else {
            // TODO: Gather pixels of slave processes
            data2pixels(data_odd, pixels);
        }
        plot(pixels);
        #endif
    }

    printf("Converge after %d iterations, elapsed time: %.6f, average computation time: %.6f\n", count-1, total_time, (double) total_time / (count-1));

    delete[] data_odd;
    //delete[] data_even;
    delete[] fire_area;

    #ifdef GUI
    delete[] pixels;
    #endif
}




int main(int argc, char *argv[]) {
    size = atoi(argv[1]);

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    rowSize = (size-2)/world_size +2;


	if (my_rank == 0) {
        #ifdef GUI
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
        glutInitWindowPosition(0, 0);
        glutInitWindowSize(window_size, window_size);
        glutCreateWindow("Heat Distribution Simulation Sequential Implementation");
        gluOrtho2D(0, resolution, 0, resolution);
        #endif
        master();
	} else {
        slave();
    }

	if (my_rank == 0){
		printf("Student ID: 117010332\n"); // replace it with your student id
		printf("Name: XU Jiale\n"); // replace it with your name
		printf("Assignment 4: Heat Distribution Simulation MPI Implementation\n");
	}

	MPI_Finalize();

	return 0;
}

