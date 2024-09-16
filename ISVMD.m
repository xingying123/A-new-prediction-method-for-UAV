function [u,u_hat,omega]= ISVMD(signal,maxAlpha,tau,tol,stopc,init_omega)

if mod(length(signal),2)>0
    signal=signal(1:end-1);
end


y = sgolayfilt(signal,8,25); 
signoise=signal-y; 

save_T = length(signal);
fs = 1/save_T;
T = save_T;
f_mir=zeros(1,T/2);
f_mir_noise=zeros(1,T/2);
f_mir(1:T/2) = signal(T/2:-1:1);
f_mir_noise(1:T/2) = signoise(T/2:-1:1);
f_mir(T/2+1:3*T/2) = signal;
f_mir_noise(T/2+1:3*T/2) = signoise;
f_mir(3*T/2+1:2*T) = signal(T:-1:T/2+1);
f_mir_noise(3*T/2+1:2*T) = signoise(T:-1:T/2+1);

f = f_mir;
fnoise=f_mir_noise;

T = length(f);
t = (1:T)/T;

udiff = tol+eps; 

omega_freqs = t-0.5-1/T;

f_hat = fftshift((fft(f)));
f_hat_onesided = f_hat;
f_hat_onesided(1:T/2) =0;
f_hat_n = fftshift((fft(fnoise)));
f_hat_n_onesided = f_hat_n;
f_hat_n_onesided(1:T/2) =0;

noisepe=norm(f_hat_n_onesided,2).^2;


N = 600;


omega_L = zeros(N, 1);

switch nargin
    case 6
        if init_omega == 0
            omega_L(1) = 0;
        else
            omega_L(1) = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,1)));
        end
    otherwise
        init_omega = 0;
        omega_L(1) = 0;
end

minAlpha=10; 
Alpha=minAlpha; 
alpha=zeros(1,1);
n = 1; 
m=0;     
SC2=0; 
l=1;  %------ the initial number of modes
bf=0;  % ----- bit flag to increase alpha
BIC=zeros(1,1);  % ------- the initial value of Bayesian index

h_hat_Temp=zeros(2, length(omega_freqs));%-initialization of filter matrix

u_hat_Temp=zeros(1,length(omega_freqs),1);%- matrix1 of modes 
u_hat_i=zeros(1, length(omega_freqs));%- matrix2 of modes

n2=0; % ---- counter for initializing omega_L


polm=zeros(2,1);  % ---- initializing Power of Last Mode index

omega_d_Temp=zeros(1,1);%-initialization of center frequencies vector1
sigerror=zeros(1,1);%initializing signal error index for stopping criteria
gamma=zeros(1,1);%----initializing gamma
normind=zeros(1,1);

while (SC2~=1)
    
    while (Alpha(1,1)<(maxAlpha+1)) 
        
        while ( udiff > tol &&  n < N ) 
            
            
                     u_hat_L(n+1,:)= (f_hat_onesided+...
                ((Alpha(1,1).^2)*(omega_freqs - omega_L(n,1)).^4).*u_hat_L(n,:)+...
                lambda(n,:)/2)./(1+(Alpha(1,1).^2)*(omega_freqs - omega_L(n,1)).^4 ...
                .*((1+(2*Alpha(1,1))*(omega_freqs - omega_L(n,1)).^2))+sum(h_hat_Temp));
            
            
            omega_L(n+1,1) = (omega_freqs(T/2+1:T)*(abs(u_hat_L(n+1, T/2+1:T)).^2)')/sum(abs(u_hat_L(n+1,T/2+1:T,1)).^2);
            
            
            lambda(n+1,:) = lambda(n,:) + tau*(f_hat_onesided...
                -(u_hat_L(n+1,:) + (((Alpha(1,1).^2)*(omega_freqs - omega_L(n,1)).^4....
                .*(f_hat_onesided - u_hat_L(n+1,:)-sum(u_hat_i)+lambda(n,:)/2)-sum(u_hat_i))...
                ./(1+(Alpha(1,1).^2)*(omega_freqs - omega_L(n,1)).^4 ))+...
                sum(u_hat_i)));
            
            
            udiff = eps;

            udiff = udiff + (1/T*(u_hat_L(n+1,:)-u_hat_L(n,:))*conj((u_hat_L(n+1,:)-u_hat_L(n,:)))')...
                / (1/T*(u_hat_L(n,:))*conj((u_hat_L(n,:)))');
            
            udiff = abs(udiff);
            
            n = n+1;
            
        end
        
            
        if abs(m-log(maxAlpha))> 1
            m=m+1;
        else
            m=m+.05;
            bf=bf+1;
        end
        if  bf>=2
            Alpha=Alpha+1;
        end
        if   Alpha(1,1)<=(maxAlpha-1)  %exp(SC1)<=(maxAlpha)
            
            if (bf ==1)
                Alpha(1,1)=maxAlpha-1;
            else
                Alpha(1,1)=exp(m);
            end
            
            omega_L=omega_L(n,1);
            
            udiff = tol+eps; % update step
            
            temp_ud = u_hat_L(n,:);%keeping the last update of obtained mode
            
            n = 1; % loop counter
            
            lambda = zeros(N, length(omega_freqs));
            u_hat_L = zeros(N, length(omega_freqs));
            u_hat_L(n,:)=temp_ud;
            
        end
    end
    
      
    omega_L=omega_L(omega_L>0);
    u_hat_Temp(1,:,l)=u_hat_L(n,:);
    omega_d_Temp(l)=omega_L(n-1,1);
    alpha(1,l)=Alpha(1,1);
    Alpha(1,1)=minAlpha;
    bf=0;
    %------------------------------initializing omega_L
    if init_omega >0
        ii=0;
        while (ii<1 && n2 < 300)
            
            omega_L = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,1)));
            
            checkp=abs(omega_d_Temp-omega_L);
            
            if (size(find(checkp<0.02),2)<=0) % it will continue if difference between previous vector of omega_d and the current random omega_plus is about 2Hz
                ii=1;
            end
            n2=n2+1;
        end
        
    else
        omega_L=0;
    end
    udiff = tol+eps; % update step

    lambda = zeros(N, length(omega_freqs));
    
    gamma(l)=1;
    
    
    h_hat_Temp(l,:)=gamma(l) ./((alpha(1,l)^2)*...
        (omega_freqs - omega_d_Temp(l)).^4);
    
    %---------keeping the last desired mode as one of the extracted modes
    u_hat_i(l,:)=u_hat_Temp(1,:,l);
    
          
    
    if nargin >=5 % checking input of the function
        
        switch stopc
            
            case 1
                %-----------------In the Presence of Noise
                if size(u_hat_i,1) == 1
                    sigerror(l)=  norm((f_hat_onesided-(u_hat_i)),2)^2;
                else
                    sigerror(l)=  norm((f_hat_onesided-sum(u_hat_i)),2)^2;
                end
                
                if ( n2 >= 300 || sigerror(l) <= round(noisepe))
                    SC2=1;
                end
            case 2
                %-----------------Exact Reconstruction
                sum_u=sum(u_hat_Temp(1,:,:),3); % -- sum of current obtained modes
                
                normind(l)=(1/T) *(norm(sum_u-f_hat_onesided).^2)...
                    ./((1/T) * norm(f_hat_onesided).^2);
                if( n2 >= 300 || normind(l) <.005 )
                    SC2=1;
                end
            case 3
                %------------------Bayesian Method
                if size(u_hat_i,1) == 1
                    sigerror(l)= norm((f_hat_onesided-(u_hat_i)),2)^2;
                    
                else
                    sigerror(l)= norm((f_hat_onesided-sum(u_hat_i)),2)^2;
                    
                end
                BIC(l)=2*T*log(sigerror(l))+(3*l)*log(2*T);
                
                if(l>1)
                    if(BIC(l)>BIC(l-1))
                        SC2=1;
                    end
                end
                
            otherwise
                %------------------Power of the Last Mode
                if (l<2)
                    polm(l)=norm((4*Alpha(1,1)*u_hat_i(l,:)./(1+2*Alpha(1,1)*...
                        (omega_freqs-omega_d_Temp(l)).^2))*u_hat_i(l,:)',2);
                    polm_temp=polm(l);
                    polm(l)=polm(l)./max(polm(l));
                else
                    polm(l)=norm((4*Alpha(1,1)*u_hat_i(l,:)./(1+2*Alpha(1,1)*...
                        (omega_freqs-omega_d_Temp(l)).^2))*u_hat_i(l,:)',2);
                    polm(l)=polm(l)./polm_temp;
                end
                
                if (l>1 && (abs(polm(l)-polm(l-1))<0.001) )
                    SC2=1;
                end
        end
    else
        %------------------Power of the Last Mode
        if (l<2)
            polm(l)=norm((4*Alpha(1,1)*u_hat_i(l,:)./(1+2*Alpha(1,1)*...
                (omega_freqs-omega_d_Temp(l)).^2))*u_hat_i(l,:)',2);
            polm_temp=polm(l);
            polm(l)=polm(l)./max(polm(l));
        else
            polm(l)=norm((4*Alpha(1,1)*u_hat_i(l,:)./(1+2*Alpha(1,1)*...
                (omega_freqs-omega_d_Temp(l)).^2))*u_hat_i(l,:)',2);
            polm(l)=polm(l)./polm_temp;
        end
        
        if (l>1 && (abs(polm(l)-polm(l-1))<tol) )
            SC2=1;
        end
    end
    
    
    
    
    
    u_hat_L = zeros(N, length(omega_freqs));
    
    
    n = 1; % ----- reset the loop counter
    
    
    l=l+1; %---(number of obtained modes)+1
    m=0;
    n2=0;
end






omega = omega_d_Temp;
L=length(omega); %------number of modes


u_hat = zeros(T, L);
u_hat((T/2+1):T,:) = squeeze(u_hat_Temp(1,(T/2+1):T,:));
u_hat((T/2+1):-1:2,:) = squeeze(conj(u_hat_Temp(1,(T/2+1):T,:)));
u_hat(1,:) = conj(u_hat(end,:));


u = zeros(L,length(t));

for l = 1:L
    u(l,:)=real(ifft(ifftshift(u_hat(:,l))));
end

[omega,indic]=sort(omega);
u=u(indic,:);
%---------- remove mirror part
u = u(:,T/4+1:3*T/4);

%--------------- recompute spectrum
clear u_hat;
u_hat=zeros(save_T,L);

for l = 1:L
    u_hat(:,l)=fftshift(fft(u(l,:)))';
end

end
