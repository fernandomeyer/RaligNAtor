################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/divsufsort/divsufsort.c \
../src/divsufsort/sssort.c \
../src/divsufsort/trsort.c \
../src/divsufsort/utils.c 

OBJS += \
./src/divsufsort/divsufsort.o \
./src/divsufsort/sssort.o \
./src/divsufsort/trsort.o \
./src/divsufsort/utils.o 

C_DEPS += \
./src/divsufsort/divsufsort.d \
./src/divsufsort/sssort.d \
./src/divsufsort/trsort.d \
./src/divsufsort/utils.d 


# Each subdirectory must supply rules for building sources it contributes
src/divsufsort/%.o: ../src/divsufsort/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


