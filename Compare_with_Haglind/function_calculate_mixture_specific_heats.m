function [MW,Cp_per_kg,Cv_per_kg,gamma] = function_calculate_mixture_specific_heats(Gas_Temperature,mole_fraction_vector,properties)
                          
%Returns Molecular Weight (MW) of mixture, Specific Heats (Cp,Cv) in kJ/kg, specific heat ratio (gamma) and
%sensible enthalpy per  (difference from 298), and entropy when supplied with Gas_Temperature in Kelvin and molar fractions
%mole fraction vector may be replaced by n kmoles of reactants or products per kmol fuel

%Unpacking mole fractions,mole_fraction_vector = [Y_N2;Y_O2;Y_AR;Y_CO2;Y_H20;Y_fuel];
Y_N2 = mole_fraction_vector(1);Y_O2 = mole_fraction_vector(2);Y_AR = mole_fraction_vector(3);Y_CO2 = mole_fraction_vector(4);Y_H2O = mole_fraction_vector(5);
%Unpacking property data
AR_data = properties.AR_data; AR_MW = properties.AR_MW;CO2_data= properties.CO2_data ;CO2_breakpoint = properties.CO2_breakpoint;CO2_MW = properties.CO2_MW;
H2O_data = properties.H2O_data;H2O_breakpoint = properties.H2O_breakpoint;H2O_MW = properties.H2O_MW;N2_data = properties.N2_data;
N2_breakpoint = properties.N2_breakpoint;N2_MW = properties.N2_MW;O2_data = properties.O2_data;O2_breakpoint = properties.O2_breakpoint;
O2_MW = properties.O2_MW;

%All correlations in terms of kJ/kmol.K
T = Gas_Temperature; %in Kelvins
t = T/1000;
%Nitrogen kJ/kmol, Cp_N2 =  28.90 - 0.1571e-2*T + 0.8081e-5*(T^2)  -  2.873e-09*(T^3);%Valid till 1800 K, from Cengel
if T<= N2_breakpoint(1)
  use_index_N2=1;%use first column of coefficients
elseif (T>N2_breakpoint(1))&&(T<=N2_breakpoint(2))
  use_index_N2=2;%use second column coefficients
elseif (T>N2_breakpoint(2))&&(T<=6000)
  use_index_N2=3;%use third column coefficients  
end
Cp_N2 =           N2_data(1:5,use_index_N2)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kj/kmol.K 

%Oxygen kJ/kmol, Cp_O2 =  25.48 +  1.520e-2*T - 0.7155e-5*(T^2)  +  1.312e-09*(T^3);%Valid till 1800 K, from Cengel
if T<= O2_breakpoint(1)
    use_index_O2=1;%use first column of coefficients         
elseif (T>O2_breakpoint(1))&&(T<=O2_breakpoint(2))
    use_index_O2=2;%use second column coefficients
elseif (T>O2_breakpoint(2))&&(T<=6000)
    use_index_O2=3;%use third column coefficients  
end
Cp_O2 =           O2_data(1:5,use_index_O2)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kJ/kmol.K    

%Carbon Dioxide kJ/kmol, Cp_CO2 = 22.26 +  5.981e-2*T -  3.501e-5*(T^2)  +  7.469e-09*(T^3);%Valid till 1800 K, from Cengel
if T<= CO2_breakpoint
    use_index_CO2=1;%use first column of coefficients          
elseif (T>CO2_breakpoint)&&(T<=6000)
    use_index_CO2=2;%use second column of coefficients 
end
Cp_CO2 =           CO2_data(1:5,use_index_CO2)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kJ/kmol.K    


%H2O kJ/kmol
if T<= H2O_breakpoint(1)
     use_index_H2O=1;%use first column of coefficients           
elseif (T>H2O_breakpoint(1))&&(T<=H2O_breakpoint(2))
    use_index_H2O=2;%use second column of coefficients 
elseif (T>H2O_breakpoint(2))&&(T<=6000)
    use_index_H2O=3;%use third column of coefficients  
end
Cp_H2O =           H2O_data(1:5,use_index_H2O)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kJ/kmol.K    

%AR
Cp_AR =           AR_data(1:5,1)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kJ/kmol.K    

%Now calculating Cp of mixture
Cp_per_mole = Y_N2*Cp_N2 + Y_O2*Cp_O2 + Y_AR*Cp_AR + Y_CO2*Cp_CO2  +  Y_H2O*Cp_H2O;%kJ/kmol.K
Cv_per_mole = Cp_per_mole-8.31446262;%kJ/kmol.K
gamma = Cp_per_mole/Cv_per_mole;
MW =          Y_N2*N2_MW + Y_O2*O2_MW + Y_AR*AR_MW + Y_CO2*CO2_MW  +  Y_H2O*H2O_MW;%kg/kmol
Cp_per_kg = Cp_per_mole/MW;%kJ/kg.K
Cv_per_kg = Cv_per_mole/MW;%kJ/kg.K

end