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

int main(int argc, char *argv[])
{
	adi_initComponents();
	pwrSetup();
	uartSetup(9600);


	uint16_t local_Id;
	uint16_t local_Net;
	uint32_t local_UniqueID;

	DEBUG_MESSAGE("starting...");
	if(LocalRead(&local_Id, &local_Net, &local_UniqueID) != MESH_OK){
		DEBUG_MESSAGE("LOCAL READ ERROR");
	}else{
		DEBUG_MESSAGE("id: %d net: %d UID: %X", local_Id, local_Net, local_UniqueID);
	}

	
	return 0;
}

