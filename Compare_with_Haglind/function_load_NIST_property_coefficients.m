function [properties] = function_load_NIST_property_coefficients


    %Loading Property Files
    load NIST_Property_Data.mat 
    properties.AR_data = AR_data;
    properties.AR_MW =  AR_MW;
    properties.AR_h_formation = AR_h_formation;
    properties.CH4_data = CH4_data;
    properties.CH4_breakpoint = CH4_breakpoint;
    properties.CH4_MW = CH4_MW ;
    properties.CH4_h_formation = CH4_h_formation;
    properties.CO2_data = CO2_data;
    properties.CO2_breakpoint = CO2_breakpoint;
    properties.CO2_MW = CO2_MW;
    properties.CO2_h_formation = CO2_h_formation;
    properties.H2_data = H2_data;
    properties.H2_breakpoint = H2_breakpoint;
    properties.H2_MW = H2_MW;
    properties.H2_h_formation = H2_h_formation;
    properties.H2O_data = H2O_data;
    properties.H2O_breakpoint = H2O_breakpoint;
    properties.H2O_MW = H2O_MW;
    properties.H2O_h_formation = H2O_h_formation;
    properties.N2_data = N2_data;
    properties.N2_breakpoint = N2_breakpoint;
    properties.N2_MW = N2_MW;
    properties.N2_h_formation = N2_h_formation;
    properties.O2_data = O2_data;
    properties.O2_breakpoint=O2_breakpoint;
    properties.O2_MW = O2_MW;
    properties.O2_h_formation = O2_h_formation;


end