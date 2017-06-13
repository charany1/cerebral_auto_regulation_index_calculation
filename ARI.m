function [a,b,c,k,v_fit,cc,ari]=ARI(v,bp,sampling_frequency,cbp,cv)

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
  keyboard
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
  v_fit(i,1:12)=v1;%v1 is current velocity wave after arma containing N_fit data points
  %v_template=zeros(9,12);
   %ee=zeros(1,9);
   
  correlation_coefficients = zeros(1,9);
  for n=1:1:9    %loop to calculate V_template(sub loop inside this loop) for all the combination of KDT
      
    K=[0.2,0.4,0.6,0.8,0.9,0.94,0.96,0.97,0.98];% gain
    D=[1.6,1.5,1.15,.9,.75,.65,.55,.52,.5];%damping factor
    T=[2,2,2,2,1.9,1.6,1.2,0.87,0.65];%time Constant
    
    
    v_template=zeros(9,12);
    %normalize current_bp_block to get dp
    keyboard
    dp = (current_bp_block - cbp)/(1-(12/cbp));
    
    %Here value of 
    for h=1:1:12% loop to calculate v_template (Refer Appendix A of the paper 
        %"Spontaneous fluctuations in cerebral blood flow regulation: contribution of PaCO2" for detail of below written lines)
        dp=((current_bp_block(h)-cbp))/(1-(12/cbp));%dp (normalized pressure change)
        keyboard
        %[x1,x2]=find_x1_x2(dp,2,T(1,n),D(1,n));% x1,x2
        x2 = get_x2_t(h);
        mV=1+dp-(K(1,n)*x2);%relative velocity
        v_template(1,h)=cv*mV;   %v_Template
    end
    %[YOGI:] check here that v_template is indeed 9x12
    keyboard
    %onle correlation coefficient for ARI 1 is coming not for entire 9 KDT
    %combinition
    correff=corrcoef(v_template(n,:),v1);% correlation coefficient b/w v_template nd V1(velocity value with N_fit)
    correlation_coeff_with_current_template=correff(1,2);% correlation coefficient from output
    correlation_coefficients(1,n)=correlation_coeff_with_current_template;%looped qoutput for 9KDT combinition
  end
  
  [max_correlation_coeff,index_of_max_correlation_coeff] = max(correlation_coefficients);
  %since ARIs are same as indexes ,i.e. {1,2,3... 9 } index for which we
  %got max value of correlation coefficient is the ARI we get my matching
  %v_fit with all 9 rows of v_template , where each corresponds to one ARI
  %value .
  
  keyboard
  cc = max_correlation_coeff;
  keyboard
  ari = index_of_max_correlation_coeff;
    
%plot(gain_factor_K);
%ylabel("Gain parameter K`");
end
