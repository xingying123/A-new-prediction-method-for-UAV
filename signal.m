clc
close all
signal = xlsread('data2.xlsx');   
plot(signal);
set(gcf,'color','w')
global tau tol stopc                  
tau = 0;                             
tol = 1e-6;                           
stopc = 4;                            
lb = 1000;                            
ub = 5000;                           
dim = 1;                              
SearchAgents_no = 30;                
Max_iter = 1000;                      
value = 1;

fobj = @(x)getObjValue(x,signal,value); 

[Best_score , Best_pos, Curve, process] = RBMO(SearchAgents_no,Max_iter,lb,ub,dim,fobj);

figure
plot(Curve);

figure
plot(process)
ylabel('MaxAlpha')
set(gcf,'color','w')

display(['SVMD最佳参数为: ', num2str(round(Best_pos))]);  
display(['SVMD最佳适应度值为 : ', num2str(Best_score)]);  

maxAlpha = Best_pos(1);                %  compactness of mode
[imf, uhat, omega] = ISVMD(signal, maxAlpha, tau, tol, stopc);
xlswrite('RBMO-SVMD分解结果.xlsx',imf);

plot_imf(imf,signal);