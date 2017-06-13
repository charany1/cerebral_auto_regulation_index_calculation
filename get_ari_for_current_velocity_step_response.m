function[ari] = get_ari_for_current_velocity_step_response(current_velocity_step_reponse,v_tempalates)

greatest_correlation_coeff_till_now = -Inf;
ari = -Inf ;

for current_ari_for_matching=1:1:9
    velocity_template_for_current_ari = v_tempalates(current_ari_for_matching,:);
    correlation_matrix = corrcoef(current_velocity_step_reponse,velocity_template_for_current_ari);
    correlation_coeff = abs(correlation_matrix(1,2));
    
    if correlation_coeff > greatest_correlation_coeff_till_now
        greatest_correlation_coeff_till_now = correlation_coeff;
        ari = current_ari_for_matching;
    end
    
end

%validation check
if ari == -Inf 
    sprintf("ARI still -Inf , something is wrong with ARI calculation");
end

end
