HH
==
% Set initial conditions 
%TIME
simulationTime = 100; %in milliseconds
dt=.01;
time=0:dt:simulationTime;

%PARAMETERS
C_m=1 %nF

g_bar_K=36; %mS/cm^2
g_bar_Na=120; %mS/sm^2
g_bar_L=0.3; %mS/cm^2

E_K = -12; %mV
E_Na=115; %mV
E_L=10.6; %mV
V_rest=-70; %mV


%COEFF and VARIABLES
V=0
alpha_n=.01*((10-V)./(exp((10-V)/10)-1)); 
beta_n=.125*exp(-V/80); 
alpha_m=.1*((25-V)./(exp((25-V)/10)-1)); 
beta_m=4*exp(-V/18); 
alpha_h=.07*exp(-V/20); 
beta_h=1./(exp((30-V)/10)+1); 

n(1)=alpha_n/(alpha_n+beta_n); 
m(1)=alpha_m/(alpha_m+beta_m); 
h(1)=alpha_h/(alpha_h+beta_h); 




for i=1:numel(time)-1; %Compute coefficients, currents, and derivates at each time step
   
    %answer to Question 1
    %I(i)=0;
    
    %Answer to Question 2.  To prevent your model neuron from dying, you decide to stimulate the
%model neuron with a step pulse of 5microAmpes for 0.5 ms!!

    %if i<=50;
        %I(i)=0.005; %mAmp
    %else
        %I(i)=0;
    %end
  
    %answer to Question 3
    %I(i)=50;
    
    %calculate the coefficients at each step
   
    alpha_n(i)=.01*((10-V(i))/(exp((10-V(i))/10)-1));
    beta_n(i)=.125*exp(-V(i)/80);
    alpha_m(i)=.1*((25-V(i))/(exp((25-V(i))/10)-1));
    beta_m(i)=4*exp(-V(i)/18);
    alpha_h(i)=.07*exp(-V(i)/20);
    beta_h(i)=1/(exp((30-V(i))/10)+1);
   
   
    %calculate the currents with different gating variables at different
    %time points
    I_Na=(m(i)^3)*g_bar_Na*h(i)*(V(i)-E_Na); 
    I_K=(n(i)^4)*g_bar_K*(V(i)-E_K); 
    I_L=g_bar_L*(V(i)-E_L); 
    I_ion=I(i)-I_K-I_Na-I_L;

   
    %calculate the derivatives using Euler approximation
    V(i+1)=V(i)+dt*I_ion/C_m;
    n(i+1)=n(i)+dt*(alpha_n(i)*(1-n(i))-beta_n(i)*n(i)); 
    m(i+1)=m(i)+dt*(alpha_m(i)*(1-m(i))-beta_m(i)*m(i)); 
    h(i+1)=h(i)+dt*(alpha_h(i)*(1-h(i))-beta_h(i)*h(i)); 

end

%Resting potential to -70mv
V=V-70; %mV


%plot Voltage
figure(1)
p1=plot(time,V,'k-');
title('Change in voltage');
legend({'Voltage'});
ylabel('Voltage (mV)');
xlabel('Time (ms)');


%plot conductances for Potassium and Sodium
figure (2)
p2=plot(time,g_bar_K*n.^4,'b');
hold on
p3=plot(time,g_bar_Na*(m.^3).*h,'r');
%p4=plot(time,g_bar_L,'k'); It seems that g_L is not important fot the
%report
legend([p2, p3], 'Conductance of Potassium', 'Conductance of Sodium');
axis([0 100 -5 40]);
ylabel('Conductance');
xlabel('Time (ms)');
title('Conductance for Potassium and Sodium ions');

%plot g_Na/g_K
figure (3)
p5=plot(time,(g_bar_Na*m.^3.*h)./(g_bar_K*n.^4));
title('Sodium conductance relative to Potassium conductance');
xlabel('Time (ms)');
ylabel('Sodium relative to Potassium conductance');

%plot m, n and h
figure(4)
p6=plot(time, m, 'g-');
title('Changes in gating variables');
hold on
p7=plot(time, n, 'b-');
p8=plot(time, h, 'r-');
legend([p6, p7, p8], 'm', 'n', 'h');
xlabel('Time (ms)');
ylabel('m, n and h');



