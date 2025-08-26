function reactant_mole_fractions = function_get_reactant_composition(intake_pressure,intake_temperature,relative_humidity)


     tmptr_array =  [0.01  5        10      13     15       20      25     30       35     40      45   50] + 273.15;%Deg Kelvin
     pg_array = [0.6117 0.8725  1.2281  1.497   1.7057 2.3392    3.1698 4.2469 5.6291 7.3851 9.5953 12.352]*1000;%Pa

     p_vapor = (relative_humidity/100).*interp1(tmptr_array,pg_array,intake_temperature);%partial pressure of water in Pascals
     p_dry_air = intake_pressure-p_vapor;
     p_N2 = 0.78084*p_dry_air;
     p_O2 = 0.20946*p_dry_air;
     p_AR = 0.00934*p_dry_air;
     p_CO2 = p_dry_air-p_N2-p_O2-p_AR;

     Y_N2 =  p_N2/intake_pressure;
     Y_O2 =  p_O2/intake_pressure;
     Y_AR=   p_AR/intake_pressure;
     Y_CO2 = p_CO2/intake_pressure;
     Y_H2O = p_vapor/intake_pressure;

     reactant_mole_fractions = [Y_N2;Y_O2;Y_AR;Y_CO2;Y_H2O];

 

end