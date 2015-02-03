################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../HLMTools/Cluster.c \
../HLMTools/HLMCopy.c \
../HLMTools/LAdapt.c \
../HLMTools/LBuild.c \
../HLMTools/LFoF.c \
../HLMTools/LGCopy.c \
../HLMTools/LGList.c \
../HLMTools/LGPrep.c \
../HLMTools/LLink.c \
../HLMTools/LMerge.c \
../HLMTools/LNewMap.c \
../HLMTools/LNorm.c \
../HLMTools/LPlex.c \
../HLMTools/LSubset.c 

OBJS += \
./HLMTools/Cluster.o \
./HLMTools/HLMCopy.o \
./HLMTools/LAdapt.o \
./HLMTools/LBuild.o \
./HLMTools/LFoF.o \
./HLMTools/LGCopy.o \
./HLMTools/LGList.o \
./HLMTools/LGPrep.o \
./HLMTools/LLink.o \
./HLMTools/LMerge.o \
./HLMTools/LNewMap.o \
./HLMTools/LNorm.o \
./HLMTools/LPlex.o \
./HLMTools/LSubset.o 

C_DEPS += \
./HLMTools/Cluster.d \
./HLMTools/HLMCopy.d \
./HLMTools/LAdapt.d \
./HLMTools/LBuild.d \
./HLMTools/LFoF.d \
./HLMTools/LGCopy.d \
./HLMTools/LGList.d \
./HLMTools/LGPrep.d \
./HLMTools/LLink.d \
./HLMTools/LMerge.d \
./HLMTools/LNewMap.d \
./HLMTools/LNorm.d \
./HLMTools/LPlex.d \
./HLMTools/LSubset.d 


# Each subdirectory must supply rules for building sources it contributes
HLMTools/%.o: ../HLMTools/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


