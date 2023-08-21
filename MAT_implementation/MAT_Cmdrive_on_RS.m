%% So Here I am trying to do the unthinkable. Trimming down the Tarnaud model to its bare essentials.

% Created: 27.July.2023 by Mithun. 
% Trying to use a capacitive drive to recreate the result of Plaksin (2014) Fig. 2b
% The Capacitance value is approximated using a sinusoidal function matching the simulated values by funPES function written by Tarnaud (2019)

tic;
%% Initialization

% Only saving the Q. 
LowMemoryMode = 1;

% Function Call parameters
Tsim = 0.085; % 0.085;
USpstart = 0.01; %0.01;
USpd = 0.05; %0.05; % For Plaksin (2016) Fig 9A
ESpstart = 0.3;
ESpd = 0.4;
ESipa = 0.0; % 0.017 is Equivalent to 0.65 nA for NEURON model 

USfreq = 690e3;
USipa = 3200; % 
aBLS = 32e-9;
fBLS = 1;

% Removing the BLS-related variables and functions. Only HH model and its solution is required. 
% RS neuron model parameters
Gna = 560;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;				% Na nernst potential (mV)
Gk = 60;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Gm = 0.75;				% Maximal conductance of the slow non-inactivating K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gl = 0.205;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -70.3;				% Leak nernst potential (mV)
VT = -56.2;             % Spike treshold adjustment parameter (mV)
taumax = 608*10^(-3);   % Decay time constant slow non-inactivating K+ (s)
Vm0 = -71.9;            % rest membrane potential (mV)
delta = 1.26*10^(-9);   % Initial gap between leaflets (charge) (m)


Cm0 = 0.01; % Resting Membrane Capacitance (F/m2) - note the units (SI). The equations are adjusted according to this. 

% Rate equations (RS neuron model)
am = @(V) -1000*0.32*((-4)*double((V-VT-13)==0)+...
    ((V-VT-13)/(double((V-VT-13)~=0)*exp(-((V-VT-13)/4))-1)));% Rate constant alpha_m [1/s]
bm = @(V) 1000*0.28*(5*double((V-VT-40)==0)+...
    (V-VT-40)/(double((V-VT-40)~=0)*exp(((V-VT-40)/5))-1)); % Rate constant beta_m [1/s]
an = @(V) -1000*0.032*((-5)*double((V-VT-15)==0)+...
    (V-VT-15)/(double((V-VT-15)~=0)*exp(-((V-VT-15)/5))-1));% Rate constant alpha_n [1/s]
bn = @(V) 1000*0.5*exp(-(V-VT-10)/40);                    % Rate constant beta_n [1/s]
ah = @(V) 1000*0.128*exp(-((V-VT-17)/18));                % Rate constant alpha_h [1/s]
bh = @(V) (1000*4)/(1+exp(-((V-VT-40)/5)));               % Rate constant beta_h [1/s]
pinf = @(V) 1/(1+exp(-((V+35)/10)));                      % Rest p-value [-]
taup = @(V) taumax/(3.3*exp((V+35)/20)+exp(-(V+35)/20));  % Time-constant [s] for p 

USstep = @(t) double(t>=USpstart&t<=USpd+USpstart); 
ESstep = @(t) double(t>=ESpstart&t<=ESpd+ESpstart); 
ESi = @ (t) ESipa*ESstep(t);  % [A/m^2]

%% Initial Conditions and Timespan

% For the neuron model
m0 = am(Vm0)/(am(Vm0)+bm(Vm0));
n0 = an(Vm0)/(an(Vm0)+bn(Vm0));
h0 = ah(Vm0)/(ah(Vm0)+bh(Vm0));
p0 = pinf(Vm0);

% Tupdate 0.025 ms, just like NEURON's default time step. 
Tupdate = 25*10^(-6);          % Update period:  every Tupdate (s) the charge Cm*Vm is updated (Plaksin et al., 2014, appendix)

% Defining more dt values for different purposes. Check when are each of them being used
dtUS = 0.025/ USfreq;
dtES = 10*10^(-6);              % ES discr. time (s)

% This is important: They moved from Vm to Q here. 
Q0 = Cm0*(10^(-3)*Vm0);      % Update charge [C]

% Initialization for HH simulation
Y0 = [Q0 m0 n0 p0 h0]';

% This marks the end of initialization. And this much runs without any error.  

%% Running simulation: or the SOLVER

% Looks like they are preparing to run the solver with a different output vector. Let's see. 
TvaluesY = [0]; 
U0 = Y0; Y=Y0';

if LowMemoryMode
    Y = Y0(1);
    Cm = [Cm0];
end

tCURRENT = 0;

% The "great" while loop begins...
while tCURRENT < Tsim

	% Setting the Zone and the time window (tNICEc) in which the simulation needs to be done. 
	if tCURRENT < USpstart
		ZONE = 1;
		tNICEc = [tCURRENT, min(USpstart,Tsim)]; % Complete (c) tNICE interval
	elseif tCURRENT >= USpstart && tCURRENT < USpstart+USpd
		ZONE = 2;
		tNICEc = [tCURRENT, min(USpstart+USpd,Tsim)];
	elseif tCURRENT >= USpstart+USpd
		ZONE = 1;
		tNICEc = [tCURRENT,Tsim];
	end

	% BLS calculation is over, now we go into the HH calculation
	if ZONE == 1
		% Just set the Cm and dt
		%CmR = @(t) Cm0;
		dt = dtES;
    else
        dt = dtUS;
    end
    CmR = @(t) FunCm_drive(Cm0, USfreq, t, ZONE); %% sinusoidal capacitance 
    tNICE=tNICEc;

	% Electrical stimulation (usually set to zero) 
	discontsES = [ESpstart,ESpstart+ESpd];
	
	% Solving HH equations using ode113
	OdeOpts=odeset('MaxStep',dt,'AbsTol',1e-3,'RelTol',1e-3);
	[t,U] = ode113(@(t,U) SimplNICERSFS(ESi,t,U(1),U(2),U(3),U(4),U(5),CmR,Gna,Vna,Gk,Vk,Gm,Gl,Vl,am,bm,an,bn,pinf,taup,ah,bh),tNICE,U0,OdeOpts); % --> External function SimplNICERSFS
	% Append the results
	TvaluesY = horzcat(TvaluesY,t(2:end)');
	% Also mark the matrix indices to which the values to be appended
	InTb = size(Y,1)+1; InTe = InTb+length(t)-2;

    Y(InTb:InTe,1) = U(2:end,1);
    if ZONE == 1
        Cm = horzcat(Cm, ones(length(t(2:end)), 1)' * Cm0);
    else
        Cm = horzcat(Cm, CmR(t(2:end)'));
    end

	% Initialization for the next cycle
	U0 = U(end,:)'; Q0 = U(end,1);

	% The following block checks if the length of the output vector matches with the length of the time vector. 
	if size(Y,1) ~= size(TvaluesY,2)
	    error('Error: size of Y-matrix and TvaluesY not equal');
	end
	tCURRENT = t(end);
    toc;
end



%% Plotting the results

% This is important: They moved from Vm to Q here. 
Potential = Y(:,1) * 1e3 ./ Cm';      % Potential evaluated from Charge and Capacitance
figure(1);
plot(TvaluesY * 1e3, Potential, 'linewidth', 2 );
hold on;
plot(TvaluesY*1e3, Cm*1000 + Vm0, 'r')
xlabel("Time (ms");
ylabel("Potential (mV)");



% Just plotting the charge for now. 
figure(2);
Charge = 10^5*Y(:,1); % Charge [nC/cm^2]
plot( TvaluesY * 1000, Charge, 'linewidth', 2 );
xlabel("Time (ms)");
ylabel("Charge (nC/cm2)");

% Save time (ms), charge (nC/cm2) and Potential (mV) into a file. 
save("../PyNEURON_implementation/mat_output.mat", "TvaluesY", "Charge", "Potential");

