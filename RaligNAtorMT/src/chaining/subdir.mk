################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/chaining/chain2dim.c \
../src/chaining/prsqualint.c \
../src/chaining/rbtree.c 

OBJS += \
./src/chaining/chain2dim.o \
./src/chaining/prsqualint.o \
./src/chaining/rbtree.o 

C_DEPS += \
./src/chaining/chain2dim.d \
./src/chaining/prsqualint.d \
./src/chaining/rbtree.d 


# Each subdirectory must supply rules for building sources it contributes
src/chaining/%.o: ../src/chaining/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF "$(@:%.o=%.d)" -MT "$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


