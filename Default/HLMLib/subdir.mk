################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../HLMLib/LCMap.c \
../HLMLib/LGBase.c \
../HLMLib/LModel.c \
../HLMLib/LPCalc.c \
../HLMLib/LPMerge.c \
../HLMLib/LUtil.c \
../HLMLib/LWMap.c 

OBJS += \
./HLMLib/LCMap.o \
./HLMLib/LGBase.o \
./HLMLib/LModel.o \
./HLMLib/LPCalc.o \
./HLMLib/LPMerge.o \
./HLMLib/LUtil.o \
./HLMLib/LWMap.o 

C_DEPS += \
./HLMLib/LCMap.d \
./HLMLib/LGBase.d \
./HLMLib/LModel.d \
./HLMLib/LPCalc.d \
./HLMLib/LPMerge.d \
./HLMLib/LUtil.d \
./HLMLib/LWMap.d 


# Each subdirectory must supply rules for building sources it contributes
HLMLib/%.o: ../HLMLib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


