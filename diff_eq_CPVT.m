function dY=diff_eq_CPVT(t,Y,con,mutation)

Vm=Y(1); %[mV] Membrane potential
qa=Y(2); %Activation gating variable of I_st
qi=Y(3); %Inactivation gating variable of I_st
dT=Y(4); %Inactivation gate of I_CaT
fT=Y(5); %Activation gate of I_CaT
pa=Y(6); %Activation gating variable of I_Kr
pi=Y(7); %Inactivation gating variable of I_Kr
xs=Y(8); %Activation gating variable of I_Ks
fL12=Y(9); %Inactivation gate for I_(C aL,1.2)
dL12=Y(10); %Activation gate for I_(CaL,1.2)
fL13=Y(11); %Inactivation gate for I_(CaL,1.3)
dL13=Y(12); %Activation gate for I_(CaL,1.3)
fCa=Y(13); %Ca2+-dependent inactivation gating variable for I_(CaL,1.2)and I_(CaL,1.2)
r=Y(14); %Activation gating variable of I_toand I_sus
q=Y(15); %Inactivation gating variable of I_to
m15=Y(16); %Activation gating variable of Nav1.5
h15=Y(17); %Fast inactivation gating variable of Nav1.5
j15=Y(18); %Slow inactivation gating variable of Nav1.5
m11=Y(19); %Activation gating variable of Nav1.1
h11=Y(20); %Fast inactivation gating variable of Nav1.1
j11=Y(21); %Slow inactivation gating variable of Nav1.1
y=Y(22); %Activation gating variable of I_f
Cai=Y(23); %[mM] Intracellular Ca2+ concentration or Ca2+ concentration in the cytosol
CajSR=Y(24); %[mM] Ca2+ concentration in the JSR
CanSR=Y(25); %[mM] Ca2+ concentration in the NSR
Casub=Y(26); %[mM] Ca2+ concentration in the subspace
f_TC=Y(27); %Fractional occupancy of the troponin Ca2+ site by [Ca2+]i
f_TMC=Y(28); %Fractional occupancy of the troponin Mg2+ site by [Ca2+]i
f_TMM=Y(29); %Fractional occupancy of the troponin Mg2+ site by Mg2+
f_CMs=Y(30); %Fractional occupancy of calmodulin by [Ca2+]sub
f_CMi=Y(31); %Fractional occupancy of calmodulin by [Ca2+]i
f_CQ=Y(32); %Fractional occupancy of calsequestrin by [Ca2+]rel
R=Y(33); %Fraction of reactivated (closed) RyR channels
OO=Y(34); %Open fraction of RyR channels
S=Y(35); %Inactive fraction of RyR channels
RI=Y(36); %Fraction of RyR inactivated channels
w=Y(37); %IKACh ACh and voltage-dependent gating variable
cAMP=Y(38); %cAMP
PLB=Y(39); %PLB phosphorylation level
A=Y(40); %The density of regulatory units with bound Ca2+ and adjacent weak cross-bridges
TT=Y(41); %The density of regulatory units with bound Ca2+ and adjacent strong cross-bridge
U=Y(42); %The density of regulatory units without bound but with adjacent strong cross-bridge
SL=Y(43); %Sarcomere length

%Mutations

if mutation==3 %hHCN4-573x
    cAMP_temp=0;
else 
    cAMP_temp=cAMP;
end

%Reversal potential
ENa = con.E_T*log(con.Nao/con.Nai);
EK=con.E_T*log(con.Ko/con.Ki);
EKs=con.E_T*log((con.Ko+0.12*con.Nao)/(con.Ki+0.12*con.Nai));
ECa=(con.E_T/2)*log(con.Cao/Casub);

%PKA acitivity4

PKA=update_PKA(cAMP);

%ATP-ADP
ATP=con.ATP_max*((con.kATP*(cAMP*100/con.cAMPb)^con.n_ATP)/(con.k_ATP05+(cAMP*100/con.cAMPb)^con.n_ATP)-con.K_ATPmin)/100;

%cAMP activity

k_ibmx = -0.8730.*con.IBMX.^0.8395./(4.0550.^0.8395+con.IBMX.^0.8395)+1; 

k_iso=0.1599*(con.ISO^1.5/(76.5441^0.6238+con.ISO^1.5)); % ISO adult - change the power (was 0.6238) 

k_CCh=0.0146*(con.CCh^1.4402/(51.7331^1.4402+con.CCh^1.4402));
k_1=con.K_ACI+con.K_AC/(1+exp((con.K_Ca-con.k_bCM*f_CMi/(con.k_fCM*(1-f_CMi)))/con.K_AC_Ca));
k_2=k_ibmx*265.3512*(cAMP^5.7343)/(24.7290^6.7343+cAMP^6.7343);
k_3=(con.k_PKA*cAMP^((con.n_PKA-1)))/(con.k_PKA_cAMP^(con.n_PKA)+cAMP^(con.n_PKA));

dcAMP = (k_iso*(ATP*0.6*1e3)+k_1*(ATP*0.6*1e3)-k_2*cAMP-k_3*cAMP-k_CCh*(ATP*0.6*1e3))/60000;

%PLB activity

k_4=((con.k_PLBp)*(con.PKA_PLB*PKA)^(con.n_PLB))/((con.k_PKA_PLB)^(con.n_PLB)+(con.PKA_PLB*PKA)^(con.n_PLB)); %***change n_PLB with ISO (20 nM)
k_5=con.k_PP1*con.PP1*PLB/(con.k_pp1_PLB+PLB); %***
dPLB=(k_4-k_5)/60000;

%4-aminopyridine-sensitive currents, I_toand I_sus
It0=con.g_to*(Vm-EK)*q*r; %nA
Isus=con.g_sus*(Vm-EK)*r; %nA
q_inf=1/(1+exp((Vm+49)/13));
r_inf=1/(1+exp(-(Vm-19.3)/15));
tau_q=(6.06+39.102/(0.57*exp(-0.08*(Vm+44))+0.065*exp(0.1*(Vm+45.93))))/0.67;
tau_r=(2.75+14.40516/(1.037*exp(0.09*(Vm+30.61))+0.369*exp(-0.12*(Vm+23.84))))/0.303;

%Ca2+ background current, IbCa
Ib_Ca=con.g_bCa*(Vm-ECa);  %nA

%K+ background current, IbK
Ib_K=con.g_bK*(Vm-EK);  %nA

%Na+ background current, IbNa
Ib_Na=con.g_bNa*(Vm-ENa);  %nA

%L-type channel current, ICaL
b_CaL=-0.2152+1.6913*PKA^10.0808/(0.8836^10.0808+PKA^10.0808); %new
ICaL12=con.g_CaL12*(1+b_CaL)*(Vm-con.E_CaL)*dL12*fL12*fCa; %new
ICaL13=con.g_CaL13*(1+b_CaL)*(Vm-con.E_CaL)*dL13*fL13*fCa; %new
ICaL=ICaL12+ICaL13;  %nA

dL12_inf=1/(1+exp(-(Vm+3)/5));
fL12_inf=1/(1+exp((Vm+36)/4.6));
dL13_inf=1/(1+exp(-(Vm+13.5)/6));
fL13_inf=1/(1+exp((Vm+35)/7.3));
fCa_inf=con.K_mfCa/(con.K_mfCa+);

%alpha_dL=-28.39*(Vm+35)/(exp(-(Vm+35)/2.5)-1)-84.9*Vm./(exp(-0.208*Vm)-1);
%beta_dL=11.43*(Vm-5)./(exp(0.4*(Vm-5))-1);
if abs(Vm)<=0.001
    alpha_dL = -28.39.*(Vm+35.0)./(exp(-(Vm+35.0)./2.5)-1.0)+408.173;
elseif abs(Vm+35)<=0.001
    alpha_dL = 70.975-84.9.*Vm./(exp(-0.208.*Vm)-1.0);
elseif abs(Vm)>0.001 && abs(Vm+35)>0.001
    alpha_dL = -28.39.*(Vm+35)./(exp(-(Vm+35)/2.5)-1)-84.9*Vm./(exp(-0.208.*Vm)-1);
end

if abs(Vm-5)<=0.001
    beta_dL = 28.575;
else
    beta_dL = 11.43*(Vm-5.0)./(exp(0.4*(Vm-5.0))-1.0);
end

tau_dL=2000/(alpha_dL+beta_dL);
tau_fL=(7.4+45.77*exp(-0.5*(Vm+28.1)*(Vm+28.1)/(11*11)));
tau_fCa=fCa_inf/con.alpha_fCa;

%T-type Ca2+ current, ICaT
ICaT=con.g_CaT*(Vm-con.E_CaT)*dT*fT;  %nA

dT_inf=1/(1+exp(-(Vm+26)/6));
fT_inf=1/(1+exp((Vm+61.7)/5.6));
tau_dT=1/(1.068*exp((Vm+26.3)/30)+1.068*exp(-(Vm+26.3)/30));
tau_fT=1/(0.0153*exp(-(Vm+61.7)/83.3)+0.015*exp((Vm+61.7)/15.38));

%Hyperpolarization-activated, funny current, If
IfNa=0.3833*con.g_If*(Vm-ENa)*y;  %nA
IfK=0.6167*con.g_If*(Vm-EK)*y;  %nA
If=IfNa+IfK;  %nA

K_if=25.3403;
K_05if=18.1115;
n_if=9.2453; %was 9.6383
V_shift=K_if*(cAMP_temp^n_if)/(K_05if^(n_if)+cAMP_temp^n_if)-18.1040;
y_inf=1/(1+exp((Vm+104.2-V_shift)/16.3)); %new

tau_y=1.5049/(exp(-(Vm+590.3)*0.01094)+exp((Vm-85.1)/17.2));

%K1 current, IK1
xk1inf = 1/(1+exp(0.070727*(Vm-EK)));
IK1=con.g_K1*xk1inf*(con.Ko/(con.Ko+0.228880))*(Vm-EK);  %nA

%Kr current, IKr
IKr=con.g_Kr*(Vm-EK)*pa*pi;  %nA
pa_inf=1/(1+exp(-(Vm+21.173694)/9.757086));
pi_inf=1/(1+exp((Vm+20.758474-4)/(19)));
tau_pa=0.699821/(0.003596*exp((Vm)/15.339290)+0.000177*exp(-(Vm)/25.868423));
tau_pi=0.2+0.9*1/(0.1*exp(Vm/54.645)+0.656*exp(Vm/106.157));

%Ks current, IKs
IKs=con.g_Ks*(Vm-EKs)*(xs^2);  %nA
xs_inf=1/(1+exp(-(Vm-20.876040)/11.852723));
tau_xs=1000/(13.097938/(1+exp(-(Vm-48.910584)/10.630272))+exp(-Vm/35.316539));

%Sodium current, INa
FNa=((9.52e-02)*exp((-6.3e-2)*(Vm+34.4))/(1+1.66*exp(-0.225*(Vm+63.7))))+8.69e-2;
hs11=(1-FNa)*h11+FNa*j11;
hs15=(1-FNa)*h15+FNa*j15;

%INa11=con.g_Na11*(m11^3)*hs11*Vm*con.Nao*con.F/(con.E_T*1000)*(exp((Vm-ENa)/con.E_T)-1)/(exp(Vm/con.E_T)-1);
%INa15=con.g_Na15*(m15^3)*hs15*Vm*con.Nao*con.F/(con.E_T*1000)*(exp((Vm-con.ENa15)/con.E_T)-1)/(exp(Vm/con.E_T)-1);


if abs(Vm)>0.005
    INa11 = con.g_Na11*m11^3*hs11*Vm*con.Nao*con.F/(con.E_T*1000)*(exp((Vm-ENa)/con.E_T)-1)/(exp(Vm/con.E_T)-1);  %nA
else
    INa11 = con.g_Na11*m11^3*hs11*Vm*con.Nao*(con.F/1000)*((exp((Vm-ENa)*1/(con.E_T))-1.0));  %nA
end

if abs(Vm)>0.005
    INa15 = con.g_Na15*m15^3*hs15*Vm*con.Nao*con.F/(con.E_T*1000)*(exp((Vm-con.ENa15)/con.E_T)-1)/(exp(Vm/con.E_T)-1);  %nA
else
    INa15 = con.g_Na15*m15^3*hs15*Vm*con.Nao*(con.F/1000)*((exp((Vm-con.ENa15)*1/(con.E_T))-1.0));  %nA
end


INa=INa11+INa15;  %nA

m11_inf=1/(1+exp(-(Vm+36.097331-5)/5))^(1/3);
h11_inf=1/(1+exp((Vm+56)/3));
j11_inf=h11_inf;
m15_inf=1/(1+exp(-(Vm+45.213705)/7.219547))^(1/3);
h15_inf=1/(1+exp(-(Vm+62.578120 )/(-6.084036)));
j15_inf=h15_inf;
tau_m11=1000*((0.6247e-03/(0.832*exp(-0.335*(Vm+56.7))+0.627*exp(0.082*(Vm+65.01))))+0.0000492);
tau_h11=1000*(((3.717e-06*exp(-0.2815*(Vm+17.11)))/(1+0.003732*exp(-0.3426*(Vm+37.76))))+0.0005977);
tau_j11=1000*(((0.00000003186*exp(-0.6219*(Vm+18.8)))/(1+0.00007189*exp(-0.6683*(Vm+34.07))))+0.003556);
tau_m15=tau_m11;
tau_h15=tau_h11;
tau_j15=tau_j11;

%Na+ - K+ pump current, INaK
INaK=con.I_NaKmax*(((con.Ko).^1.2)./((con.K_mK).^1.2+(con.Ko).^1.2))*((con.Nai.^1.3)./...
    ((con.K_mNa).^1.3+(con.Nai).^1.3))./(1.0+exp(-(Vm-ENa+120)/30));  %nA
  
%Na+-Ca2+ exchanger current, INCX

d0=1+(con.Cao/con.K_co)*(1+exp(con.Q_co*Vm/con.E_T))+(con.Nao/con.K_1no)*(1+(con.Nao/con.K_2no)*(1+con.Nao/con.K_3no)); 
k43=con.Nai/(con.K_3ni+con.Nai);
k41=exp(-con.Q_n*Vm/(2*con.E_T));
k34=con.Nao/(con.K_3no+con.Nao);
k21=(con.Cao/con.K_co)*exp(con.Q_co*Vm/con.E_T)/d0;
k23=(con.Nao/con.K_1no)*(con.Nao/con.K_2no)*(1+con.Nao./con.K_3no)*exp(-con.Q_n*Vm./(2*con.E_T))./d0;
k32=exp(con.Q_n*Vm./(2*con.E_T));
x1=k34*k41*(k23+k21)+k21*k32*(k43+k41);
di=1+(Casub/con.K_ci)*(1+exp(-con.Q_ci*Vm/con.E_T)+con.Nai/con.K_cni)+(con.Nai/con.K_1ni)*...
    (1+(con.Nai/con.K_2ni)*(1+con.Nai/con.K_3ni));
k12= (Casub/con.K_ci)*exp(-con.Q_ci*Vm/con.E_T)/di;
k14=(con.Nai/con.K_1ni)*(con.Nai/con.K_2ni)*(1+con.Nai/con.K_3ni)*exp(con.Q_n*Vm./(2*con.E_T))./di;
x2=k43*k32*(k14+k12)+k41*k12*(k34+k32);
x3=k43*k14*(k23+k21)+k12*k23*(k43+k41);
x4=k34*k23*(k14+k12)+k21*k14*(k34+k32);

INaCa=con.K_NaCa*(k21*x2-k12*x1)/(x1+x2+x3+x4);  %nA

%Sustained inward current, Ist
I_st=con.g_st*(Vm-con.E_st)*qa*qi;  %nA
qa_inf=1/(1+exp(-(Vm+67)/5));
alpha_qa=1/(0.15*exp(-Vm/11)+0.2*exp(-Vm/700));
beta_qa=1/(16*exp(Vm/8)+15*exp(Vm/50));
tau_qa=1/(alpha_qa+beta_qa);
alpha_qi=0.15*1/(3100*exp((Vm+10)/13)+700.3*exp((Vm+10)/70));
beta_qi=0.15*1/(95.7*exp(-(Vm+10)/10)+50*exp(-(Vm+10)/700))+0.000229/(1+exp(-(Vm+10)/5));
qi_inf=alpha_qi/(alpha_qi+beta_qi);
tau_qi=1/(alpha_qi+beta_qi);

%Acetylcholine-activated K+ current, IKACh
I_KACh=con.C*con.g_KACh_max*(Vm-EK)*w;  %nA
beta_w=0.001*12.32/(1+0.0042/(con.CCh*10^(-6)));
alpha_w=0.001*17*exp(0.0133*(Vm+40));
w_inf=beta_w/(alpha_w+beta_w);
tau_w=1/(alpha_w+beta_w);
a_w = w_inf/tau_w;
b_w = (1-w_inf)/tau_w;
dw = a_w*(1-w)-b_w*w;

%Ca2+ fluxes in the SR

%Ryanodine Receptor

jSRCarel=con.k_s*OO*(CajSR-Casub);
kCaSR=con.MaxSR-(con.MaxSR-con.MinSR)/(1+(con.EC_50SR/CajSR)^con.HSR);
k_oCa = con.k_oCa_max*(con.RyR_min-con.RyR_max*PKA^con.n_RyR/(con.k_05Ry^con.n_RyR+PKA^con.n_RyR)+1);
koSRCa=k_oCa/kCaSR;
kiSRCa=con.k_iCa*kCaSR;
dR=(con.k_im*RI-kiSRCa*Casub*R)-(koSRCa*(Casub)^con.n*R-con.k_om*OO);
dOO=(koSRCa*(Casub)^con.n*R-con.k_om*OO)-(kiSRCa*Casub*OO-con.k_im*S);
dS=(kiSRCa*Casub*OO-con.k_im*S)-(con.k_om*S-koSRCa*(Casub)^con.n*RI);
dRI=(con.k_om*S-koSRCa*(Casub)^con.n*RI)-(con.k_im*RI-kiSRCa*Casub*R);

%Intracellular Ca2+ flux
fPLB=2.9102*PLB^9.5517/(0.2763^9.5517+PLB^9.5517)+0.4998; %+0.4998 was missing
j_up=(con.P_up*fPLB)/(1+con.K_up/Cai);


j_Cadiff=(Casub-Cai)/con.tho_difCa; %Ca2+ diffusion flux from submembrane space to myoplasm, jCa,dif
j_tr=(CanSR-CajSR)/con.tho_tr; %Ca2+ flux between network and junctional SR compartments, jtr


    % = Force
    Ve = 0; % 
    dSL = -Ve;
    NXB = (SL-con.SL_lo)/2*1e3*(TT+U)*con.Nc;
    K_Ca = con.FK0+con.Fkl*(NXB^con.FN)/(con.FK05^con.FN+NXB^con.FN);
    k_l = con.Fkl/K_Ca;
    dA = con.Fkl*Cai*(1-A-TT-U)-(con.Ff+k_l)*A+(con.Fg0+con.Fg1*Ve)*TT;
    dTT = con.Ff*A-(con.Fg0+con.Fg1*Ve+k_l)*TT+con.Fkl*Cai*U;
    dU = k_l*TT-(con.Fg0+con.Fg1*Ve+con.Fkl*Cai)*U;    


%gating variables (For any gating variable g with steady state g)
dqa=(qa_inf-qa)/tau_qa;
dqi=(qi_inf-qi)/tau_qi;
ddT=(dT_inf-dT)/tau_dT;
dfT=(fT_inf-fT)/tau_fT;
dpa=(pa_inf-pa)/tau_pa;
dpi=(pi_inf-pi)/tau_pi;
dxs=(xs_inf-xs)/tau_xs;
dfL12=(fL12_inf-fL12)/tau_fL;
ddL12=(dL12_inf-dL12)/tau_dL;
dfL13=(fL13_inf-fL13)/tau_fL;
ddL13=(dL13_inf-dL13)/tau_dL;
dfCa=(fCa_inf-fCa)/tau_fCa;
dr=(r_inf-r)/tau_r;
dq=(q_inf-q)/tau_q;
dm15=(m15_inf-m15)/tau_m15;
dh15=(h15_inf-h15)/tau_h15;
dj15=(j15_inf-j15)/tau_j15;
dm11=(m11_inf-m11)/tau_m11;
dh11=(h11_inf-h11)/tau_h11;
dj11=(j11_inf-j11)/tau_j11;
dy=(y_inf-y)/tau_y;
%dw=(w_inf-w)/tau_w;

%Ca2+ buffering
%df_TC=con.k_fTC*Cai*(1-f_TC)-con.k_bTC*f_TC;
df_TC = con.Fkl*Cai*(1-A-TT)-k_l*(A+TT); % if using force expression
df_TMC=con.k_fTMC*Cai*(1-f_TMC-f_TMM)-con.k_bTMC*f_TMC;
df_TMM=con.k_fTMM*con.Mgi*(1-f_TMC-f_TMM )-con.k_bTMM*f_TMM;
df_CMi=con.k_fCM*Cai*(1-f_CMi)-con.k_bCM*f_CMi;
df_CMs=con.k_fCM*Casub*(1-f_CMs)-con.k_bCM*f_CMs;
df_CQ=con.k_fCQ*CajSR*(1-f_CQ)-con.k_bCQ*f_CQ;   

%Dynamics of Ca2+ concentrations in cell compartments
dCai=((j_Cadiff*con.V_sub-j_up*con.V_nSR)/con.V_i)-(con.CM_tot*(df_CMi)+con.TC_tot*(df_TC)+con.TMC_tot*(df_TMC));
dCasub=(-(ICaL+ICaT+Ib_Ca-2*INaCa)/(2*con.F/1000)+jSRCarel*con.V_jSR)/con.V_sub-j_Cadiff-con.CM_tot*df_CMs;
dCajSR=j_tr-jSRCarel-con.CQ_tot*df_CQ;
dCanSR=j_up-j_tr*con.V_jSR/con.V_nSR;


%Current
%I_KACh=0;
I=INa+ICaL+ICaT+IKr+IKs+If+It0+Isus+IK1+INaK+INaCa+Ib_Na+Ib_K+Ib_Ca+I_KACh;

%voltage [mV]

dVm=-(I+I_st)/con.C; %[nA/nF]=[mV/ms]

%Derivatives

dY=zeros(43,1);

dY(1)=dVm;
dY(2)=dqa;
dY(3)=dqi;
dY(4)=ddT;
dY(5)=dfT;
dY(6)=dpa;
dY(7)=dpi;
dY(8)=dxs;
dY(9)=dfL12;
dY(10)=ddL12;
dY(11)=dfL13;
dY(12)=ddL13;
dY(13)=dfCa;
dY(14)=dr;
dY(15)=dq;
dY(16)=dm15;
dY(17)=dh15;
dY(18)=dj15;
dY(19)=dm11;
dY(20)=dh11;
dY(21)=dj11;
dY(22)=dy;
dY(23)=dCai;
dY(24)=dCajSR;
dY(25)=dCanSR;
dY(26)=dCasub;
dY(27)=df_TC;
dY(28)=df_TMC;
dY(29)=df_TMM;
dY(30)=df_CMs;
dY(31)=df_CMi;
dY(32)=df_CQ;
dY(33)=dR;
dY(34)=dOO;
dY(35)=dS;
dY(36)=dRI;
dY(37)=dw;
dY(38)=dcAMP;
dY(39)=dPLB;
dY(40)=dA;
dY(41)=dTT;
dY(42)=dU;
dY(43)=dSL;

disp(t/1000);