function [a,b,c,k,v_fit,cc]=ARI(v,bp,sampling_frequency,cbp,cv)

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
cc=zeros((number_of_moving_windows),9);
for i=1:1:(number_of_moving_windows)
  lp=119+i; %last point of the moving block
  %keyboard
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
  v1=zeros(1,12);
  for p=1:1:12 % here , Nfit was in seconds so , we have to map that to number of data points accordingly
      vel=(current_bp_block(p,1)*a(i)) + (current_diff_block_with_zeros_prepended((p+1),1)*b(i)) +( current_diff_block_with_zeros_prepended(p,1)*c(i));
      v1(p)=vel;
  end
  v_fit(i,1:12)=v1;
  %v_template=zeros(9,12);
   %ee=zeros(1,9);
  for n=1:1:9
      
    K=[0.2 0.4 0.6 0.8 0.9 0.94 0.96 0.97 0.98];
    D=[1.6 1.5 1.15 .9 .75 .65 .55 .52 .5];
    T=[2 2 2 2 1.9 1.6 1.2 0.87 0.65];
    %v_template=zeros(9,12);
    for h=1:1:12
        dp=current_bp_block(h)-cbp;
        [x1,x2]=find_x1_x2(dp,2,T(1,n),D(1,n));
         mV=1+dp-(K(1,n)*x2(1,1));
         v_template(1,h)=cv*mV;   
    end
    correff=corrcoef(v_template,v1);
    dd=correff(1,2);
    ee(n)=dd;
  end
cc(i,1:9)=ee;
        plot(gain_factor_K);
ylabel("Gain parameter K`");
end
