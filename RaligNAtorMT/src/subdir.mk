################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/align.c \
../src/alignOnline.c \
../src/alignSLink.c \
../src/alignThread.c \
../src/alignThreadInit.c \
../src/alphabet.c \
../src/init.c \
../src/resultsProcessing.c \
../src/search.c \
../src/streamHandling.c 

OBJS += \
./src/align.o \
./src/alignOnline.o \
./src/alignSLink.o \
./src/alignThread.o \
./src/alignThreadInit.o \
./src/alphabet.o \
./src/init.o \
./src/resultsProcessing.o \
./src/search.o \
./src/streamHandling.o 

C_DEPS += \
./src/align.d \
./src/alignOnline.d \
./src/alignSLink.d \
./src/alignThread.d \
./src/alignThreadInit.d \
./src/alphabet.d \
./src/init.d \
./src/resultsProcessing.d \
./src/search.d \
./src/streamHandling.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF "$(@:%.o=%.d)" -MT "$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


