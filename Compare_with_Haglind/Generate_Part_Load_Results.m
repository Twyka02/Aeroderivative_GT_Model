clearvars
close all

%Started 7/29/2024 to generate part-load results to compare with Haglind
%Need to change trb map files within function_split_shaft_model

trb_filename =  'Scaled_2Stage_HighPqPTrb_4p00_kg_s'; %TwoStgTurbine_HighPqPTurbine but 4 kg/s bleed
parasitic_losses = 4.0;%kg/s

error_threshold = 1.5;%Percent, had to increase it since RPM=3600 not chosen even at 80% load 
design_fuel_rate  = 1.934;%kg/s

N = 10;%10% to 100% load
power_vector = zeros(N,1);
torque_vector = zeros(N,1);
RPM_vector = zeros(N,1);
efficiency_vector = zeros(N,1);
mdot_vector = zeros(N,1);
PR_vector = zeros(N,1);
exh_tmptr_vector = zeros(N,1);
trg_in_tmptr_vector = zeros(N,1);
exh_p_vector = zeros(N,1);

for index = 2:10
    fuel_flow_rate = design_fuel_rate*(index/10);
    [operating_power,operating_torque, effective_RPM,operating_efficiency,operating_mdot,operating_PR,...
        operating_exhaust_tmptr_deg_C,operating_trb_inlet_tmptr_deg_C,operating_exhaust_P_kPa] = function_split_shaft_model(fuel_flow_rate,trb_filename,error_threshold,parasitic_losses);
    power_vector(index,1)=      operating_power;
    torque_vector(index,1)=     operating_torque;
    RPM_vector(index,1)=        effective_RPM;
    efficiency_vector(index,1)= operating_efficiency;
    mdot_vector(index,1)=       operating_mdot;
    PR_vector(index,1)=         operating_PR;
    exh_tmptr_vector(index,1)=  operating_exhaust_tmptr_deg_C;
    trg_in_tmptr_vector(index,1)= operating_trb_inlet_tmptr_deg_C;
    exh_p_vector(index,1)=      operating_exhaust_P_kPa;

end

PL = power_vector/max(power_vector);
eta_PL = (0.975*PL)./( 0.975*PL + (1-0.975)*( (1-0.43) +0.43*PL.*PL )  );
generator_power = power_vector.*eta_PL;
generator_load = generator_power/max(generator_power);

save Results_2Stage_HighPqPTrb_4p00_kg_s.mat