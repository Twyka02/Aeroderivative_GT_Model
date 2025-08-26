function [enthalpy_per_kmol_fuel] = function_calculate_mixture_enthalpy(Gas_Temperature,mole_fraction_vector,properties,fuel_type,fuel_temperature) 
                          
%Returns sensible enthapn given mole fraction vector (or mole vector or molar rate vector), fuel
%type and gas and fuel temperature

%Unpacking mole fractions,mole_fraction_vector = [kmol_N2;kmol_O2;kmol_AR;kmol_CO2;Y_H20;Y_fuel];
kmol_N2 = mole_fraction_vector(1);kmol_O2 = mole_fraction_vector(2);kmol_AR = mole_fraction_vector(3);kmol_CO2 = mole_fraction_vector(4);kmol_H2O = mole_fraction_vector(5);
%Unpacking property data
AR_data = properties.AR_data; AR_MW = properties.AR_MW;AR_h_formation = properties.AR_h_formation;CH4_data = properties.CH4_data;CH4_breakpoint = properties.CH4_breakpoint;
CH4_MW  = properties.CH4_MW ;CH4_h_formation = properties.CH4_h_formation ;CO2_data= properties.CO2_data ;CO2_breakpoint = properties.CO2_breakpoint;CO2_MW = properties.CO2_MW;
CO2_h_formation = properties.CO2_h_formation;H2_data = properties.H2_data;H2_breakpoint = properties.H2_breakpoint;H2_MW = properties.H2_MW;H2_h_formation = properties.H2_h_formation;
H2O_data = properties.H2O_data;H2O_breakpoint = properties.H2O_breakpoint;H2O_MW = properties.H2O_MW;H2O_h_formation = properties.H2O_h_formation;N2_data = properties.N2_data;
N2_breakpoint = properties.N2_breakpoint;N2_MW = properties.N2_MW;N2_h_formation = properties.N2_h_formation;O2_data = properties.O2_data;O2_breakpoint = properties.O2_breakpoint;
O2_MW = properties.O2_MW;O2_h_formation = properties.O2_h_formation;

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
Cp_N2 =           N2_data(1:5,use_index_N2)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kJ/kmol.K 
delta_ho_N2 =      1000*N2_data(:,use_index_N2)'*[t (t^2)/2 (t^3)/3 (t^4)/4 -1/t 1 0 -1]'; %A*t + B*t2/2 + C*t3/3 + D*t4/4 − E/t + F − H, note G is zero, units are kJ/kmol  
s_N2 =            N2_data(:,use_index_N2)'*[log(t) t (t^2)/2 (t^3)/3 -0.5/(t^2) 0 1 0]';%S° = A*ln(t) + B*t + C*t2/2 + D*t3/3 − E/(2*t2) + G, note F and H are zero, units j/mol.K

%Oxygen kJ/kmol, Cp_O2 =  25.48 +  1.520e-2*T - 0.7155e-5*(T^2)  +  1.312e-09*(T^3);%Valid till 1800 K, from Cengel
if T<= O2_breakpoint(1)
    use_index_O2=1;%use first column of coefficients         
elseif (T>O2_breakpoint(1))&&(T<=O2_breakpoint(2))
    use_index_O2=2;%use second column coefficients
elseif (T>O2_breakpoint(2))&&(T<=6000)
    use_index_O2=3;%use third column coefficients  
end
Cp_O2 =           O2_data(1:5,use_index_O2)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kJ/kmol.K    
delta_ho_O2 =      1000*O2_data(:,use_index_O2)'*[t (t^2)/2 (t^3)/3 (t^4)/4 -1/t 1 0 -1]'; %A*t + B*t2/2 + C*t3/3 + D*t4/4 − E/t + F − H, note G is zero, units are kJ/kmol  
s_O2 =            O2_data(:,use_index_O2)'*[log(t) t (t^2)/2 (t^3)/3 -0.5/(t^2) 0 1 0]';%S° = A*ln(t) + B*t + C*t2/2 + D*t3/3 − E/(2*t2) + G, note F and H are zero, units j/mol.K

%Carbon Dioxide kJ/kmol, Cp_CO2 = 22.26 +  5.981e-2*T -  3.501e-5*(T^2)  +  7.469e-09*(T^3);%Valid till 1800 K, from Cengel
if T<= CO2_breakpoint
    use_index_CO2=1;%use first column of coefficients          
elseif (T>CO2_breakpoint)&&(T<=6000)
    use_index_CO2=2;%use second column of coefficients 
end
Cp_CO2 =           CO2_data(1:5,use_index_CO2)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kJ/kmol.K    
delta_ho_CO2 =      1000*CO2_data(:,use_index_CO2)'*[t (t^2)/2 (t^3)/3 (t^4)/4 -1/t 1 0 -1]'; %A*t + B*t2/2 + C*t3/3 + D*t4/4 − E/t + F − H, note G is zero, units are kJ/kmol  
s_CO2 =            CO2_data(:,use_index_CO2)'*[log(t) t (t^2)/2 (t^3)/3 -0.5/(t^2) 0 1 0]';%S° = A*ln(t) + B*t + C*t2/2 + D*t3/3 − E/(2*t2) + G, note F and H are zero, units j/mol.K

%H2O kJ/kmol
if T<= H2O_breakpoint(1)
     use_index_H2O=1;%use first column of coefficients           
elseif (T>H2O_breakpoint(1))&&(T<=H2O_breakpoint(2))
    use_index_H2O=2;%use second column of coefficients 
elseif (T>H2O_breakpoint(2))&&(T<=6000)
    use_index_H2O=3;%use third column of coefficients  
    
end
Cp_H2O =           H2O_data(1:5,use_index_H2O)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kJ/kmol.K    
delta_ho_H2O =      1000*H2O_data(:,use_index_H2O)'*[t (t^2)/2 (t^3)/3 (t^4)/4 -1/t 1 0 -1]'; %A*t + B*t2/2 + C*t3/3 + D*t4/4 − E/t + F − H, note G is zero, units are kJ/kmol  
s_H2O =            H2O_data(:,use_index_H2O)'*[log(t) t (t^2)/2 (t^3)/3 -0.5/(t^2) 0 1 0]';%S° = A*ln(t) + B*t + C*t2/2 + D*t3/3 − E/(2*t2) + G, note F and H are zero, units j/mol.K

%AR
Cp_AR =           AR_data(1:5,1)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kJ/kmol.K    
delta_ho_AR =      1000*AR_data(:,1)'*[t (t^2)/2 (t^3)/3 (t^4)/4 -1/t 1 0 -1]'; %A*t + B*t2/2 + C*t3/3 + D*t4/4 − E/t + F − H, note G is zero, units are kJ/kmol  
s_AR =            AR_data(:,1)'*[log(t) t (t^2)/2 (t^3)/3 -0.5/(t^2) 0 1 0]';%S° = A*ln(t) + B*t + C*t2/2 + D*t3/3 − E/(2*t2) + G, note F and H are zero, units j/mol.K

%%           Fuels     %%
%CH4 kJ/kmol
t = fuel_temperature/1000;
if fuel_temperature<= CH4_breakpoint
    use_index_CH4=1;          
elseif (fuel_temperature>CH4_breakpoint)&&(fuel_temperature<=6000)
    use_index_CH4=2; 
end
Cp_CH4 =            CH4_data(1:5,use_index_CH4)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kJ/kmol.K    
delta_ho_CH4 =      1000*CH4_data(:,use_index_CH4)'*[t (t^2)/2 (t^3)/3 (t^4)/4 -1/t 1 0 -1]'; %A*t + B*t2/2 + C*t3/3 + D*t4/4 − E/t + F − H, note G is zero, units are kJ/kmol  
s_CH4 =             CH4_data(:,use_index_CH4)'*[log(t) t (t^2)/2 (t^3)/3 -0.5/(t^2) 0 1 0]';%S° = A*ln(t) + B*t + C*t2/2 + D*t3/3 − E/(2*t2) + G, note F and H are zero, units j/mol.K

%Hydrogen kJ/kmol
if fuel_temperature<= H2_breakpoint(1)
    use_index_H2=1;          
elseif (fuel_temperature>H2_breakpoint(1))&&(fuel_temperature<=H2_breakpoint(2))
    use_index_H2=2; 
elseif (fuel_temperature>H2_breakpoint(2))&&(fuel_temperature<=6000)
    use_index_H2=3;
end
Cp_H2 =           H2_data(1:5,use_index_H2)'*[1 t t^2 t^3 1/t^2]'; %A + B*t + C*t2 + D*t3 + E/t2, units are kJ/kmol.K    
delta_ho_H2 =      1000*H2_data(:,use_index_H2)'*[t (t^2)/2 (t^3)/3 (t^4)/4 -1/t 1 0 -1]'; %A*t + B*t2/2 + C*t3/3 + D*t4/4 − E/t + F − H, note G is zero, units are kJ/kmol  
s_H2 =            H2_data(:,use_index_H2)'*[log(t) t (t^2)/2 (t^3)/3 -0.5/(t^2) 0 1 0]';%S° = A*ln(t) + B*t + C*t2/2 + D*t3/3 − E/(2*t2) + G, note F and H are zero, units j/mol.K

%defaults for no fuel, fuel_type==0
h_formation_fuel = 0;
delta_ho_fuel = 0;

if fuel_type==1%pure methane
    Cp_fuel = Cp_CH4 ; %Valid till 1500 K
    MW_fuel = CH4_MW;%pure methane
    delta_ho_fuel = delta_ho_CH4;
    s_fuel =  s_CH4;
    h_formation_fuel = CH4_h_formation; % Note NIST formation data is in MJ/kmol
elseif fuel_type==2%pure hydrogen
    Cp_fuel =  Cp_H2;%Valid till 1800 K
    MW_fuel =  H2_MW;%pure Hydrogen
    delta_ho_fuel = delta_ho_H2;
    s_fuel =  s_H2;
    h_formation_fuel = 0;
elseif fuel_type==3%F-76
    Cp_fuel = Cp_CH4 ; %Valid till 1500 K
    MW_fuel = CH4_MW;%pure methane
    delta_ho_fuel = delta_ho_CH4;
    s_fuel =  s_CH4;
    h_formation_fuel = -514.0790; % Note NIST formation data is in MJ/kmol
end

%Calculating enthalpy and entropy
sensible_enthalpy_per_mole = kmol_N2*delta_ho_N2 + kmol_O2*delta_ho_O2 + kmol_AR*delta_ho_AR + kmol_CO2*delta_ho_CO2  +  kmol_H2O*delta_ho_H2O + 1*delta_ho_fuel;%kJ/kmol.K, only sensible enthalpy
chemical_enthalpy_per_mole = (kmol_N2*0 + kmol_O2*0 + kmol_AR*0 + kmol_CO2*CO2_h_formation  +  kmol_H2O*H2O_h_formation + 1*h_formation_fuel)*1000;%kJ/kmol.K, only chemical enthalpy,% Note NIST formation data is in MJ/kmol
enthalpy_per_kmol_fuel = sensible_enthalpy_per_mole + chemical_enthalpy_per_mole;
entropy_per_mole =  kmol_N2*s_N2 + kmol_O2*s_O2 + kmol_AR*s_AR + kmol_CO2*s_CO2  +  kmol_H2O*s_H2O;%J/mol.K


end