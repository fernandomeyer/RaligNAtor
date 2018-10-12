################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/afAlphabet.c \
../src/afCons.c \
../src/afStreamHandling.c \
../src/sufconstruct.c 

OBJS += \
./src/afAlphabet.o \
./src/afCons.o \
./src/afStreamHandling.o \
./src/sufconstruct.o 

C_DEPS += \
./src/afAlphabet.d \
./src/afCons.d \
./src/afStreamHandling.d \
./src/sufconstruct.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


