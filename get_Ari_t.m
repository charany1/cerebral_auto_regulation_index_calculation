function [a,b,c,k,v_fit,cc,ari]=get_Ari_t(v,bp,sampling_frequency,critical_blood_pressure,critical_velocity)

%we want to derive ARIs at every 0.5 sec i.e. frequency = (1/0.5) i.e. 2Hz
downsampling_factor = sampling_frequency/2;

downsampled_velocity=downsample(v,downsampling_factor);
%keyboard
downsampled_bp=downsample(bp,downsampling_factor);
%keyboard
data_segemnt_size=length(downsampled_velocity);
%keyboard
diff_of_bp_and_velocity=(downsampled_bp-downsampled_velocity);%difference between bp and velocity
%keyboard
number_of_moving_windows = data_segemnt_size - 119;
a=zeros((number_of_moving_windows),1);
b=zeros((number_of_moving_windows),1);
c=zeros((number_of_moving_windows),1); 
%keyboard
gain_factor_K = zeros((number_of_moving_windows),1);

%N_fit represents number of data points for template matching with Tieck's
%standard model
N_fit = 12;

% Contains first N_fit data points from velocity step response for template
% matching for number_of_moving_windows
v_fit=zeros((number_of_moving_windows),N_fit);
ari_values = zeros(1,number_of_moving_windows);
for i=1:1:(number_of_moving_windows)
  lp=119+i; %last point of the moving block
%  keyboard
  current_bp_block=downsampled_bp(i:lp);%moving block
  %current_dp=((mean(current_bp_block(1:12)))-cbp)
  %keyboard
  current_diff_block=diff_of_bp_and_velocity(i:lp);
  %keyboard
  
  %we added two zeros in beginning of current_diff_block so that ,
  %when n =1, current_diff_block(n-1)=0   
  %and when n=1&2,current_diff_block(n-2)=0 
  current_diff_block_with_zeros_prepended=[0,0,current_diff_block']';%difference moving block
  %keyboard
  %%expected dimesion of X 300x4
  X=[ones(120,1) current_bp_block current_diff_block_with_zeros_prepended(2:121) current_diff_block_with_zeros_prepended(1:120)];
  %keyboard
  %%expected dimesion of  300x1
  vb=downsampled_velocity(i:lp);
  %keyboard
  %% to find cofficients for this data segemnt , 
  %z represents coeffcient matrix
  % expected dim of z = 4x1
  z=regress(vb,X);
  %intercept_weight = z(1);
  a(i) = z(2);
  b(i) = z(3);
  c(i) = z(4);
  %keyboard
  km=(1-a(i))/(1+b(i)+c(i));
  gain_factor_K(i) = km;
  k= gain_factor_K;
  %keyboard
  %%finding N-fit sized sample of velocity step response as per ARMA
  %%co-efficients a(i),b(i),c(i)
  sampled_velocity_step_response=zeros(1,12);
  for p=1:1:12 % here , Nfit was in seconds so , we have to map that to number of data points accordingly
      vel=(current_bp_block(p,1)*a(i)) + (current_diff_block_with_zeros_prepended((p+1),1)*b(i)) +( current_diff_block_with_zeros_prepended(p,1)*c(i));
      sampled_velocity_step_response(1,p)=vel;
  end
  %keyboard
  %tranposing current_bp_block since get_v_templates expects it as row
  %vector
  current_bp_block = current_bp_block'
  
 v_templates = get_v_templates(current_bp_block,80,45);%check inside get_v_templates for values of K,D,T,f
 ari = get_ari_for_current_velocity_step_response(sampled_velocity_step_response,v_templates);
 ari_values(1,i) = ari;
 
 
end
plot(ari_values');
ylabel("Time varying ari");
end