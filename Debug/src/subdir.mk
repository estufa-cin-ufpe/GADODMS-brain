################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/LoRaMESH.c \
../src/Master.c \
../src/pwr.c \
../src/timer0.c \
../src/uart.c 

SRC_OBJS += \
./src/LoRaMESH.o \
./src/Master.o \
./src/pwr.o \
./src/timer0.o \
./src/uart.o 

C_DEPS += \
./src/LoRaMESH.d \
./src/Master.d \
./src/pwr.d \
./src/timer0.d \
./src/uart.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: CrossCore GCC ARM Embedded C Compiler'
	arm-none-eabi-gcc -g -gdwarf-2 -ffunction-sections -fdata-sections -DCORE0 -D_DEBUG -D_RTE_ -D__ADUCM3029__ -D__SILICON_REVISION__=0x100 @includes-f1055e5c1c2853d069c0dde258c52165.txt -Wall -c -mcpu=cortex-m3 -mthumb -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


