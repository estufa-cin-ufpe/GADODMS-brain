/*****************************************************************************
 * Master.c
 *****************************************************************************/

#include <sys/platform.h>
#include "adi_initialize.h"
#include "Master.h"
#include "pwr.h"
#include "timer0.h"
#include "uart.h"
#include "LoRaMESH.h"
#include <string.h>
#include "common.h"
#include <stdio.h>
#include <math.h>

float rssi_matrix[3][2];

void build_rssi(uint8_t *buffer){
	float temp;
	int i = 0;
	int j = 0;
	int k = 0;
	while(i < 24){
		temp = (float)(buffer[i] | (buffer[i+1] << 8) | (buffer[i+2] << 16) | (buffer[i+3] << 24))/100;

		i = i + 4;

		rssi_matrix[j][k++] = temp;
		if(k == 2){
			j++;
			k = 0;
		}
	}
}

void build_vector(uint8_t *buffer, double* position){
	for(int k=0;k<2;k++){
	  for(int i=1;i<5; i++){
		buffer[k*4 + i] = ((uint32_t)(position[k]*100)>>((i-1)*8))&0xFF;
	  }
	}
}

void nelder_mead_optimization(double* d, double* solution){
    //float *custo;
    double base[3][2];
    base[0][0] = 0;
    base[0][1] = 1;
    base[1][0] = 2;
    base[1][1] = 1;
    base[2][0] = 1;
    base[2][1] = 0;


    float a, b, c, D, e, f;
    a = 2*base[1][0] - 2*base[0][0];
    b = 2*base[1][1] - 2*base[0][1];
    c = d[0]*d[0] - d[1]*d[1] - base[0][0]*base[0][0] + base[1][0]*base[1][0] - base[0][1]*base[0][1] + base[1][1]*base[1][1];
    D = 2*base[2][0] - 2*base[1][0];
    e = 2*base[2][1] - 2*base[1][1];
    f = d[1]*d[1] - d[2]*d[2] - base[1][0]*base[1][0] + base[2][0]*base[2][0] - base[1][1]*base[1][1] + base[2][1]*base[2][1];
    solution[0] = (c*e - f*b) / (e*a - b*D);
    solution[1] = (c*D - a*f) / (b*D - a*e);

}

int main(int argc, char *argv[])
{
	adi_initComponents();
	pwrSetup();
	uartSetup(9600);
	uint8_t buffer_payload[24];
	uint8_t position_payload[9];
	uint8_t payload_size;
	uint16_t remote_id;


	/* EKF_p kalman_rssi[3];
	    for(int i=0;i<3;i++){
	        kalman_rssi[i].inicializar();
	        kalman_rssi[i].modelar();
	    }*/

	int i, id;
	//double dados[3][2];
	/*dados[0][0] = 19;
	dados[0][1] = 0;
	dados[1][0] = 34;
	dados[1][1] = 1;
	dados[2][0] = 19;
	dados[2][1] = 2;*/

	double d[3], n[3] = {4.29, 4.29, 4.29}, rssi_base[3] = {25.52, 25.52, 25.52}; //definir valores para n
	double ponto[2];
	double rssi;

	for(int i = 0; i < 4; i++){
		if(ReceivePacketTransp(&remote_id, buffer_payload, &payload_size, 5000) == MESH_OK){
			DEBUG_MESSAGE("id da msg: %d\n", remote_id);
			build_rssi(buffer_payload);
			for(int j = 0; j < 3; j++){
				for(int k = 0; k < 2; k++){
					DEBUG_MESSAGE("%f ", rssi_matrix[j][k]);
				}
			}
		}
	}
	for(i=0;i<3;i++){
		id = (int)(rssi_matrix[i][1]);
		rssi = rssi_matrix[i][0];
		d[i]=pow(10, ((rssi-rssi_base[id])/(10*n[id])));
	}
	nelder_mead_optimization(d, ponto);
	DEBUG_MESSAGE("Res: %f %f \n", ponto[0], ponto[1]);
	build_vector(position_payload, ponto);
	//position_payload[8] = 10;
	for(int i=0; i<9; i++){
		DEBUG_MESSAGE("%d   ", position_payload[i]);
	}
	PrepareFrameTransp(0, position_payload, sizeof(position_payload));
	if(SendPacket() == MESH_OK) DEBUG_MESSAGE("Mandei\n");
	return 0;
}

