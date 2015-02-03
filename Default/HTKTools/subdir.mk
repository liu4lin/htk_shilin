################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../HTKTools/HBuild.c \
../HTKTools/HCompV.c \
../HTKTools/HCopy.c \
../HTKTools/HDMan.c \
../HTKTools/HERest.c \
../HTKTools/HHEd.c \
../HTKTools/HInit.c \
../HTKTools/HLEd.c \
../HTKTools/HLRescore.c \
../HTKTools/HLStats.c \
../HTKTools/HList.c \
../HTKTools/HMMIRest.c \
../HTKTools/HParse.c \
../HTKTools/HQuant.c \
../HTKTools/HRest.c \
../HTKTools/HResults.c \
../HTKTools/HSGen.c \
../HTKTools/HSLab.c \
../HTKTools/HSmooth.c \
../HTKTools/HVite.c \
../HTKTools/dspr.c \
../HTKTools/dtpmc.c \
../HTKTools/tpmc.c 

OBJS += \
./HTKTools/HBuild.o \
./HTKTools/HCompV.o \
./HTKTools/HCopy.o \
./HTKTools/HDMan.o \
./HTKTools/HERest.o \
./HTKTools/HHEd.o \
./HTKTools/HInit.o \
./HTKTools/HLEd.o \
./HTKTools/HLRescore.o \
./HTKTools/HLStats.o \
./HTKTools/HList.o \
./HTKTools/HMMIRest.o \
./HTKTools/HParse.o \
./HTKTools/HQuant.o \
./HTKTools/HRest.o \
./HTKTools/HResults.o \
./HTKTools/HSGen.o \
./HTKTools/HSLab.o \
./HTKTools/HSmooth.o \
./HTKTools/HVite.o \
./HTKTools/dspr.o \
./HTKTools/dtpmc.o \
./HTKTools/tpmc.o 

C_DEPS += \
./HTKTools/HBuild.d \
./HTKTools/HCompV.d \
./HTKTools/HCopy.d \
./HTKTools/HDMan.d \
./HTKTools/HERest.d \
./HTKTools/HHEd.d \
./HTKTools/HInit.d \
./HTKTools/HLEd.d \
./HTKTools/HLRescore.d \
./HTKTools/HLStats.d \
./HTKTools/HList.d \
./HTKTools/HMMIRest.d \
./HTKTools/HParse.d \
./HTKTools/HQuant.d \
./HTKTools/HRest.d \
./HTKTools/HResults.d \
./HTKTools/HSGen.d \
./HTKTools/HSLab.d \
./HTKTools/HSmooth.d \
./HTKTools/HVite.d \
./HTKTools/dspr.d \
./HTKTools/dtpmc.d \
./HTKTools/tpmc.d 


# Each subdirectory must supply rules for building sources it contributes
HTKTools/%.o: ../HTKTools/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


