function [] = main()
 K=[0.2,0.4,0.6,0.8,0.9,0.94,0.96,0.97,0.98];% gain
 D=[1.6,1.5,1.15,.9,.75,.65,.55,.52,.5];%damping factor
 T=[2,2,2,2,1.9,1.6,1.2,0.87,0.65];%time Constant
 f = 100;
 critical_blood_pressure = 80;
 current_bp_block = evalin('base','current_bp_block');
 critical_velocity = 45;
 v_tempalates = get_v_templates(current_bp_block,critical_blood_pressure,critical_velocity,f,K,D,T);
 keyboard
 
 ari = get_ari_for_current_velocity_step_response(current_velocity_step_reponse,v_tempalates);

end
