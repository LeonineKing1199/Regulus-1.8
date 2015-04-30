################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../bucket_hash.cu \
../main.cu \
../mesh.cu \
../triangulate.cu 

CU_DEPS += \
./bucket_hash.d \
./main.d \
./mesh.d \
./triangulate.d 

OBJS += \
./bucket_hash.o \
./main.o \
./mesh.o \
./triangulate.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-7.0/bin/nvcc -O3 --use_fast_math -std=c++11 -gencode arch=compute_50,code=sm_50  -odir "." -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-7.0/bin/nvcc -O3 --use_fast_math -std=c++11 --compile --relocatable-device-code=false -gencode arch=compute_50,code=compute_50 -gencode arch=compute_50,code=sm_50  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


