function [HPC_Power_MW,T2] = function_HPC_compression(P1,P2,T1,HPC_eta,mdot_cmp,intake_air_mole_fraction,properties)

            eta_mech = 0.99;
            %Compression Process
            P_compression = [P1:1:P2]';%Compressing in steps of 1 Pa
            T_compression_isentropic = T1*ones(numel(P_compression),1);
            %compressing with variable specific heat ratio gamma
            N = numel(P_compression);
            for index = 2:N
                %get properties first
                [MW_air,Cp,Cv,gamma] = function_calculate_mixture_specific_heats(T_compression_isentropic(index-1),intake_air_mole_fraction,properties); 
                pressure_ratio =   P_compression(index)/P_compression(index-1);%for the next step increasing pressure by one Pascal
                T_compression_isentropic(index) = T_compression_isentropic(index-1) * pressure_ratio^((gamma-1)/gamma);  
            end
            T2s =  T_compression_isentropic(end);
            T2 = (T2s-T1)/HPC_eta + T1;

          %Compressor power with chemistry
          fuel_temperature = 298;%K even though there is no fuel here
          reactant_kmole_rate = intake_air_mole_fraction*(mdot_cmp/MW_air);
          h1 = function_calculate_mixture_enthalpy(T1,reactant_kmole_rate ,properties,0,fuel_temperature);%kJ/s
          h2 = function_calculate_mixture_enthalpy(T2,reactant_kmole_rate ,properties,0,fuel_temperature);%kJ/s
          HPC_Power_MW = ((h2-h1)/1000)/ eta_mech;%(kJ/kmol)*(kg/s)*(kmol/kg) = kW/1000 = MW


end