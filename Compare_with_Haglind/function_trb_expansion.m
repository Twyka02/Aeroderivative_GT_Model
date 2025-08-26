function [trb_power_MW ,T_out] = function_trb_expansion(P_in,P_out,T_in,product_kmoles_per_kmol_fuel,trb_eta,fuel_flow_rate,MW_fuel,mdot_trb)

%function to calculate expansion across a turbine given the pressure BCs
%and inlet tmptr. Uses variable specific heats and NIST chemistry

            %constants
            eta_mechanical = 0.99;
            Cp_trb = 1.16;
            fuel_temperature=298;%K

            %T_out_isen = T_in/(PR^((Cp_trb-1)/Cp_trb)); %%instead use below:
            %%%%%expansion with variable specific heat ratio gamma%%
            properties = function_load_NIST_property_coefficients;
            products_mole_fractions = product_kmoles_per_kmol_fuel/sum(product_kmoles_per_kmol_fuel);
            P_expansion = [P_in:-1:P_out]';%Compressing in steps of 10 kPa
            T_expansion_isentropic = T_in*ones(numel(P_expansion),1);
            for index = 2:numel(P_expansion)
                [MW_prod,Cp_trb_prod,Cv_prod,gamma_prod] = function_calculate_mixture_specific_heats(T_expansion_isentropic(index-1),products_mole_fractions,properties); 
                pressure_ratio =   P_expansion(index)/P_expansion(index-1);%for the next step decreasing pressure by one Pascal
                T_expansion_isentropic(index) = T_expansion_isentropic(index-1) * pressure_ratio^((gamma_prod-1)/gamma_prod);  
            end
            T_out_isen  =  T_expansion_isentropic(end);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            T_out = T_in - (T_in-T_out_isen )*trb_eta;

            %Without Chemistry
            Power_MW_cold_air_analysis = mdot_trb*Cp_trb*(T_in-T_out)/1000;
            %With Chemistry
            h3 = function_calculate_mixture_enthalpy(T_in,product_kmoles_per_kmol_fuel ,properties,0,fuel_temperature);
            h4 = function_calculate_mixture_enthalpy(T_out,product_kmoles_per_kmol_fuel ,properties,0,fuel_temperature);
            trb_power_MW=  eta_mechanical * (h3-h4)*(fuel_flow_rate/MW_fuel)/1000;%(kJ/kmol)*(kg/s)*(kmol/kg) = kW/1000 = MW


end