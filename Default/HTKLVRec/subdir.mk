################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../HTKLVRec/HDecode.mod.o \
../HTKLVRec/HDecode.orig.o \
../HTKLVRec/HLVLM.mod.o \
../HTKLVRec/HLVLM.orig.o \
../HTKLVRec/HLVModel.mod.o \
../HTKLVRec/HLVModel.orig.o \
../HTKLVRec/HLVNet.mod.o \
../HTKLVRec/HLVNet.orig.o \
../HTKLVRec/HLVRec.mod.o \
../HTKLVRec/HLVRec.orig.o 

C_SRCS += \
../HTKLVRec/HDecode.c \
../HTKLVRec/HDecode.mod.c \
../HTKLVRec/HLVLM.c \
../HTKLVRec/HLVModel.c \
../HTKLVRec/HLVNet.c \
../HTKLVRec/HLVRec-GC.c \
../HTKLVRec/HLVRec-LM.c \
../HTKLVRec/HLVRec-misc.c \
../HTKLVRec/HLVRec-outP.c \
../HTKLVRec/HLVRec-propagate.c \
../HTKLVRec/HLVRec-traceback.c \
../HTKLVRec/HLVRec.c 

OBJS += \
./HTKLVRec/HDecode.o \
./HTKLVRec/HDecode.mod.o \
./HTKLVRec/HLVLM.o \
./HTKLVRec/HLVModel.o \
./HTKLVRec/HLVNet.o \
./HTKLVRec/HLVRec-GC.o \
./HTKLVRec/HLVRec-LM.o \
./HTKLVRec/HLVRec-misc.o \
./HTKLVRec/HLVRec-outP.o \
./HTKLVRec/HLVRec-propagate.o \
./HTKLVRec/HLVRec-traceback.o \
./HTKLVRec/HLVRec.o 

C_DEPS += \
./HTKLVRec/HDecode.d \
./HTKLVRec/HDecode.mod.d \
./HTKLVRec/HLVLM.d \
./HTKLVRec/HLVModel.d \
./HTKLVRec/HLVNet.d \
./HTKLVRec/HLVRec-GC.d \
./HTKLVRec/HLVRec-LM.d \
./HTKLVRec/HLVRec-misc.d \
./HTKLVRec/HLVRec-outP.d \
./HTKLVRec/HLVRec-propagate.d \
./HTKLVRec/HLVRec-traceback.d \
./HTKLVRec/HLVRec.d 


# Each subdirectory must supply rules for building sources it contributes
HTKLVRec/%.o: ../HTKLVRec/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


