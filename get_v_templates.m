%% returns : v_templates :- matrix(9x12) where each row contains velocity templates for ari =1 ,2...9 and respective values of K,D,T
% params : 
% 1. current_bp_block : bp data for current moving window
% 2. cbp : critical_blood_pressure[VINIT , please correct if it's something else]
% 3. f : frequency 
% 4. K : gain factor (1x9)
% 5. D : Damping factor (1x9)
% 6. T : Time Constants (1x9)

%to fingure our should we pass current_bp_block or entire bp data?
function [v_templates]=get_v_templates(current_bp_block,critical_blood_pressure,critical_velocity)

 K=[0.2,0.4,0.6,0.8,0.9,0.94,0.96,0.97,0.98];% gain
 D=[1.6,1.5,1.15,.9,.75,.65,.55,.52,.5];%damping factor
 T=[2,2,2,2,1.9,1.6,1.2,0.87,0.65];%time Constant
 f = 100;

%dp : normalzed change in bp
normalized_bp_chnage = (current_bp_block - critical_blood_pressure)/(1-(12/critical_blood_pressure));


v_templates = zeros(9,12);
for ari = 1:1:9
    
    T_for_current_ari = T(1,ari)
    D_for_current_ari = D(1,ari);
    K_for_current_ari = K(1,ari);
    
    denominator=f*T_for_current_ari;
    
    
    x1_previous = 0;
    x2_previous = 0;
    dp_previous = 0;
    for t=1:1:12
        
         if t > 1
            dp_previous = normalized_bp_chnage(1,t-1);
        end
        
        %calculating x1 for current t
        x1 =    x1_previous +  ((dp_previous-x2_previous)/denominator);
        %calculating x2 for current t
        x2 = x2_previous + (((x1_previous-(2*D_for_current_ari*x2_previous)))/denominator);
        
        x1_previous = x1;
        x2_previous = x2;
       
        dp_current = normalized_bp_chnage(1,t);
        
        relative_velocity = 1+dp_current-(K_for_current_ari * x2);%relative velocity
        velocity = critical_velocity * relative_velocity;
        v_templates(ari,t) = velocity;
        
       

    end
end
end

    
