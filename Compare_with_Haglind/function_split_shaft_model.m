function [operating_power,operating_torque, effective_RPM,operating_efficiency,operating_mdot,operating_PR,...
    operating_exhaust_tmptr_deg_C,operating_trb_inlet_tmptr_deg_C,operating_exhaust_P_kPa] = function_split_shaft_model(fuel_flow_rate,trb_filename,error_threshold,parasitic_losses)


%% Compressor Maps
load Scaled_HighPqPCompr.mat scaled_corr_mdot_HPC scaled_PR_HPC scaled_eta_HPC corrected_RPM_HPC


%% Turbine Maps
eval([ 'load ' trb_filename ' scaled_PR_LPT scaled_mdot_LPT scaled_eta_LPT SF_PR_LPT SF_mdot_LPT corrected_RPM_LPT scaled_PR_HPT scaled_mdot_HPT  scaled_eta_HPT SF_PR_HPT SF_mdot_HPT corrected_RPM_HPT']);


%% End of turbine maps

%User Defined Inputs
mode = 2;%==1 for constant torque and 2 for constant speed
operating_RPM = 3600;
percent_load = 20;
operating_torque = interp1([10 100],[6 85.52],percent_load); %kN-m
T1 = 15+273.15;%Kelvin, ISO used by Haglind et.al.
T3_design = 1522.2; %Kelvin
corr_spd_design = 3600/sqrt(1110.2/T1);%LPT, based on T4
cmp_design_RPM = 8000;%From Haglind Part 2
%error_threshold = 1.0;%Percent error to tolerate when reading maps
error_exh_p = 0.5;
fuel_type = 3;%f-76
lhv = 42.798; %f-76
%fuel_flow_rate = (percent_load/100)*1.934;%design point at 100% load kg/s, from Haglind et.al.
p_ambient = 101.325;%kPa, ambient pressure
P1 = p_ambient*(1-0.01);% 1% relative inlet pressure loss from Haglind et.al.
exhaust_pressure = p_ambient*(1+0.029);% 2.9% relative inlet pressure loss from Haglind et.al.
%For comparison
lpt_power_MW = 31.207/0.975;%from Haglind et.al.
exhaust_tmptr = 533.8 + 273.15;%from Haglind et.al.
intake_RH = 0;%percent relative humidity
fuel_temperature = 298;%Deg Kelvin, Note NIST correlations are not valid below 298 deg Kelven


%% CHEMISTRY
%intake properties for reference
MW_fuel = 14.4*12.0107 + 26.9*1.00794;%F76 fuel
properties = function_load_NIST_property_coefficients;
intake_air_mole_fraction = function_get_reactant_composition(p_ambient*1000,T1,intake_RH);%Note, function requires pressure in Pa not kPa, Format = [Y_N2;Y_O2;Y_AR;Y_CO2;Y_H2O];
[MW_air,Cp_cold_air_air,Cv_air,gamma_air] = function_calculate_mixture_specific_heats(T1,intake_air_mole_fraction,properties); 

%% Start of Triple Loop
[M1, N1]=size(scaled_PR_HPC);
[M2, N2]=size(scaled_PR_HPT);
[M3, N3]=size(scaled_PR_LPT);
properties = function_load_NIST_property_coefficients;
HPC_store = zeros(M3,4);HPT_store = zeros(M3,4);LPT_store = zeros(M3,4);Power_store = zeros(M3,12);%Each LPT speed line gets a valid point
error_store = 1000*ones(M3,5);%potential solution for every speed line on the trb map

for i_HPC = 1:M1 %each i_HPC is a speed line
%for i_HPC = 61:61   %DESIGN POINT for Scaled_TwoStgTurbine_TwoStgTurbine.mat 

    HPC_RPM = corrected_RPM_HPC(i_HPC);
    for j_HPC = 1:N1
          
          %% Compression
          HPC_PR = scaled_PR_HPC(i_HPC,j_HPC);
          P2 = P1*HPC_PR;
          HPC_corrected_mdot = scaled_corr_mdot_HPC(i_HPC,j_HPC);
          mdot_cmp = HPC_corrected_mdot*(P1/P1)/sqrt(T1/T1);
          HPC_eta = scaled_eta_HPC(i_HPC,j_HPC);
          [HPC_power_MW,T2] = function_HPC_compression(P1,P2,T1,HPC_eta,mdot_cmp,intake_air_mole_fraction,properties);

          %% Combustion
          %%%%%%%%%%%% mdot_trb = 0.99*mdot_cmp+ fuel_flow_rate;% 1% bleed air, Haglind
          mdot_trb = mdot_cmp+fuel_flow_rate-parasitic_losses;% bleed air
          %%mdot_trb = 0.89*mdot_cmp+ fuel_flow_rate;%Experiment with prefix 'Exp'

          P3 = 0.97*P2;%from Haglind, 3% pressure loss
          %Combustion Chemistry
          AFR = mdot_trb/fuel_flow_rate;
          kmol_air_per_kmol_fuel = (MW_fuel*AFR)/MW_air;%numerator is mass of air per kmol fuel
          reactant_kmoles_per_kmol_fuel = kmol_air_per_kmol_fuel*intake_air_mole_fraction;%mole vector, excludes fuel
          x_r = reactant_kmoles_per_kmol_fuel;
          if fuel_type==3 %for F76
            product_kmoles_per_kmol_fuel =  [x_r(1);  x_r(2)-21.525;  x_r(3);  x_r(4)+14.8;  x_r(5)+13.45];%%%%[kmol_N2;kmol_O2;kmol_AR;kmol_CO2;kmol_H2O]; 
          end
          %Using optimization to find adiabatic product tmptr
          [reactant_enthalpy_per_mole] = 0.99*function_calculate_mixture_enthalpy(T2,reactant_kmoles_per_kmol_fuel ,properties,fuel_type,fuel_temperature);%99 percent cmbstn efficiency
           fun = @(x)(reactant_enthalpy_per_mole - function_calculate_mixture_enthalpy(x,product_kmoles_per_kmol_fuel ,properties,0,fuel_temperature))^2;
          options = optimoptions('fmincon','Display','off');
           [T_product_adiabatic,residual]  = fmincon(fun,T2+500,[],[],[],[],T2,6000,[],options);%
          if sqrt(residual)>=1
              disp('Calculation for Adiabatic Product Tmptr did not converge to less than 1 Deg C)')
              pause
          end
          T3 = T_product_adiabatic;
          %% HPT Expansion
          HPT_corrected_mdot_required = mdot_trb*sqrt(T3/T1)/(P3/P1);
          [~,i_HPT] = min(    (  (HPC_RPM*sqrt(T3_design/T3)) - corrected_RPM_HPT ).^2      );%i_HPT is the index of the speed line:
          HPT_RPM = corrected_RPM_HPT(i_HPT);
          for j_HPT = 1:N2 %Examine all HPT points at HPC RPM
                HPT_corrected_mdot = scaled_mdot_HPT (i_HPT,j_HPT);
                error_HPT_mdot = 100*abs((HPT_corrected_mdot_required-HPT_corrected_mdot)/HPT_corrected_mdot_required);
                if error_HPT_mdot<error_threshold
                     HPT_PR =  scaled_PR_HPT(i_HPT,j_HPT);
                     P4 = P3/HPT_PR;
                     HPT_eta= scaled_eta_HPT(i_HPT,j_HPT);
                     [HPT_power_MW ,T4] = function_trb_expansion(P3,P4,T3,product_kmoles_per_kmol_fuel,...
                        HPT_eta,fuel_flow_rate,MW_fuel,mdot_trb);
                     error_HPT_work = 100*(abs(HPC_power_MW-HPT_power_MW)/HPC_power_MW);
                     if error_HPT_work<error_threshold%trb work matches cmp work
                         %% LPT Expansion
                         LPT_corrected_mdot_required = mdot_trb*sqrt(T4/T1)/(P4/P1);
                         for i_LPT = 1:M3 %look at each speed line of the LPT
                             LPT_RPM = corrected_RPM_LPT(i_LPT);
                             min_error = 1000;
                             for j_LPT = 1:N3                         
                                 LPT_corrected_mdot = scaled_mdot_LPT(i_LPT,j_LPT);
                                     error_LPT_mdot = 100*abs((LPT_corrected_mdot_required-LPT_corrected_mdot)/LPT_corrected_mdot_required);
                                     if error_LPT_mdot< error_threshold
                                                    LPT_PR =  scaled_PR_LPT(i_LPT,j_LPT);
                                                    P5 = P4/LPT_PR;
                                                    LPT_eta= scaled_eta_LPT(i_LPT,j_LPT);
                                                    error_backpressure = 100*abs((P5-exhaust_pressure)/exhaust_pressure) ;
                                                    if error_backpressure  < error_threshold %valid points below
                                                        [LPT_Power_MW_shaft ,T5] = function_trb_expansion(P4,P5,T4,product_kmoles_per_kmol_fuel,...
                                                           LPT_eta,fuel_flow_rate,MW_fuel,mdot_trb);     
                                                           cycle_efficiency = 100*LPT_Power_MW_shaft/(fuel_flow_rate*lhv);
                                                           this_error = error_HPT_mdot+error_HPT_work+error_LPT_mdot+error_backpressure; 
                                                           
                                                           if this_error < error_store(i_LPT,1)
                                                               error_store(i_LPT,1) = this_error;
                                                               error_store(i_LPT,2:5) = [error_HPT_mdot error_HPT_work error_LPT_mdot error_backpressure];
                                                               HPC_store(i_LPT,:) = [ HPC_corrected_mdot HPC_PR HPC_eta HPC_RPM ];
                                                               HPT_store(i_LPT,:) = [ HPT_corrected_mdot HPT_PR HPT_eta HPT_RPM];
                                                               LPT_store(i_LPT,:) = [ LPT_corrected_mdot LPT_PR LPT_eta LPT_RPM ];
                                                               Power_store(i_LPT,:) = [ LPT_Power_MW_shaft cycle_efficiency T1 T2 T3 T4 T5 P1 P2 P3 P4 P5];
                                                           end


                                                     end %of speed_flag=0
                                     end % if LPT mdot matches
                                    
                              end%of j_LPT
                         end% of i_LPT
                           
                     end
                end
          end %end of j_HPT
        
    end%end of j_HPC
end%end of i_HPC
%Plotting

remove_index = find(error_store(:,1)==1000);
HPC_store(remove_index,:)=[];
HPT_store(remove_index,:)=[];
LPT_store(remove_index,:)=[];
Power_store(remove_index,:)=[];
error_store(remove_index,:)=[];

%CHOOSING TRB speed based on torque (RPM ranges from 800 to 360 torque = 85.52 kN-m
RPM_candidates =    LPT_store(:,4).*corr_spd_design.*sqrt(Power_store(:,6)/T1);% f.corr_spd.sqrt(To/T1)
torque_candidates = 1000*Power_store(:,1)./(2*pi*RPM_candidates/60);%in kN-m
if mode==1%constant torque mode, CHOOSING TRB speed based on torque (RPM ranges from 800 to 360 torque = 85.52 kN-m
    [closest_torque,chosen_index] = min(abs(operating_torque-torque_candidates));
elseif mode==2%constant RPM mode, RPM mode, CHOOSING TRB speed based on 3600 RPM
    [closest_RPM,chosen_index] = min(abs(operating_RPM-RPM_candidates));
end
    
operating_mdot = HPC_store(chosen_index,1)
operating_PR = HPC_store(chosen_index,2)
effective_RPM =  RPM_candidates(chosen_index)
operating_torque = torque_candidates(chosen_index)
operating_power =  Power_store(chosen_index,1)
operating_efficiency =       Power_store(chosen_index,2)
operating_trb_inlet_tmptr_deg_C =  Power_store(chosen_index,5)-273.15
operating_exhaust_tmptr_deg_C =  Power_store(chosen_index,7)-273.15
operating_exhaust_P_kPa =  Power_store(chosen_index,12)

figure(4),plot(torque_candidates, Power_store(:,2),'.',torque_candidates(chosen_index), Power_store(chosen_index,2),'x','MarkerSize',8);
xlabel('Torque [kN-m]');ylabel('Efficiency [%]'); legend('Possible Operating Points','Operating Point');
figure(5),plot(RPM_candidates, Power_store(:,2),'.',RPM_candidates(chosen_index), Power_store(chosen_index,2),'x','MarkerSize',8);
xlabel('RPM');ylabel('Efficiency [%]');legend('Possible Operating Points','Operating Point');
toc

figure(1),plot(scaled_corr_mdot_HPC',scaled_PR_HPC','.',HPC_store(chosen_index,1),HPC_store(chosen_index,2),'*');
figure(1),hold('on');xlabel('Corrected Mass Flow [kg/s]');ylabel('Pressure Ratio [None]');title('Compressor Map');
[c,h] = contour(scaled_corr_mdot_HPC,scaled_PR_HPC,scaled_eta_HPC,50);
figure(1),hold('off');

figure(2),plot(scaled_PR_HPT',scaled_mdot_HPT',':',HPT_store(chosen_index,2),HPT_store(chosen_index,1),'o','MarkerSize',8);
xlabel('Pressure Ratio [None]');ylabel('Corrected Mass Flow [kg/s]');legend('HPT Map','Operating Point');

figure(3),plot(scaled_PR_LPT',scaled_mdot_LPT',':',LPT_store(chosen_index,2),LPT_store(chosen_index,1),'o','MarkerSize',8);
xlabel('Pressure Ratio [None]');ylabel('Corrected Mass Flow [kg/s]');legend('LPT Map','Operating Point');


end