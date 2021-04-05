
TARGET = green
F90_SOURCE = ./source/linalg.f90 \
             ./source/Green.f90 \
             ./source/misc.f90 \
             ./source/main.f90 

#$(wildcard ./source/*.f90)
OBJ = $(subst .f90,.o,$(subst source,objects,$(F90_SOURCE)))
FC = gfortran
FC_FLAGS = -O3 \
           -ffast-math \
           -llapack \
           -lblas \
           -fopenmp \

MOD = mods           
           
RM = rm -rf           

all : display $(TARGET)
	

display:
	@clear
	@mkdir -p objects
	@mkdir -p mods
	
	@echo ".------------------------------------."
	@echo "| Compiling: Green's functions       |#"
	@echo "| ---------                          |#"
	@echo "|                                    |#"
	@echo "| By: Prof. Carlo R. da Cunha        |#"
	@echo "|                                    |#"
	@echo "| Created: Mar/2020                  |#"
	@echo "| Revision: Apr/2021                 |#"
	@echo "'------------------------------------'#"
	@echo "  #####################################"
	@echo ""


$(TARGET): $(OBJ)
	@echo "# Linking $@..."
	$(FC) $^ $(FC_FLAGS) -o $@ -I $(MOD)
	@echo ""

./objects/%.o: ./source/%.f90
	@echo "# Building target: $<"
	$(FC) $< $(FC_FLAGS) -c -o $@ -J $(MOD) 
	@echo ""
		
	
clean:
	@$(RM) ./objects/*.o $(TARGET) *~
	@$(RM) ./mods/*.mod $(TARGET) *~
	@rmdir objects
	@rmdir mods
	
.PHONY: all clean	
	
	
