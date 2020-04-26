function I = CPVT_equations(Y,con,mutation)

Vm=Y(:,1); %[mV] Membrane potential
qa=Y(:,2); %Activation gating variable of I_st
qi=Y(:,3); %Inactivation gating variable of I_st
dT=Y(:,4); %Inactivation gate of I_CaT
fT=Y(:,5); %Activation gate of I_CaT
pa=Y(:,6); %Activation gating variable of I_Kr
pi=Y(:,7); %Inactivation gating variable of I_Kr
xs=Y(:,8); %Activation gating variable of I_Ks
fL12=Y(:,9); %Inactivation gate for I_(CaL,1.2)
dL12=Y(:,10); %Activation gate for I_(CaL,1.2)
fL13=Y(:,11); %Inactivation gate for I_(CaL,1.3)
dL13=Y(:,12); %Activation gate for I_(CaL,1.3)
fCa=Y(:,13); %Ca2+-dependent inactivation gating variable for I_(CaL,1.2)and I_(CaL,1.2)
r=Y(:,14); %Activation gating variable of I_toand I_sus
q=Y(:,15); %Inactivation gating variable of I_to
m15=Y(:,16); %Activation gating variable of Nav1.5
h15=Y(:,17); %Fast inactivation gating variable of Nav1.5
j15=Y(:,18); %Slow inactivation gating variable of Nav1.5
m11=Y(:,19); %Activation gating variable of Nav1.1
h11=Y(:,20); %Fast inactivation gating variable of Nav1.1
j11=Y(:,21); %Slow inactivation gating variable of Nav1.1
y=Y(:,22); %Activation gating variable of I_f
Cai=Y(:,23); %[mM] Intracellular Ca2+ concentration or Ca2+ concentration in the cytosol
CajSR=Y(:,24); %[mM] Ca2+ concentration in the JSR
CanSR=Y(:,25); %[mM] Ca2+ concentration in the NSR
Casub=Y(:,26); %[mM] Ca2+ concentration in the subspace
f_TC=Y(:,27); %Fractional occupancy of the troponin Ca2+ site by [Ca2+]i
f_TMC=Y(:,28); %Fractional occupancy of the troponin Mg2+ site by [Ca2+]i
f_TMM=Y(:,29); %Fractional occupancy of the troponin Mg2+ site by Mg2+
f_CMs=Y(:,30); %Fractional occupancy of calmodulin by [Ca2+]sub
f_CMi=Y(:,31); %Fractional occupancy of calmodulin by [Ca2+]i
f_CQ=Y(:,32); %Fractional occupancy of calsequestrin by [Ca2+]rel
R=Y(:,33); %Fraction of reactivated (closed) RyR channels
OO=Y(:,34); %Open fraction of RyR channels
S=Y(:,35); %Inactive fraction of RyR channels
RI=Y(:,36); %Fraction of RyR inactivated channels
w=Y(:,37); %IKACh ACh and voltage-dependent gating variable
cAMP=Y(:,38); %cAMP
PLB=Y(:,39); %PLB phosphorylation level
A=Y(:,40); %The density of regulatory units with bound Ca2+ and adjacent weak cross-bridges
TT=Y(:,41); %The density of regulatory units with bound Ca2+ and adjacent strong cross-bridge
U=Y(:,42); %The density of regulatory units without bound but with adjacent strong cross-bridge
SL=Y(:,43); %Sarcomere length

%%

%Reversal potential
ENa = con.E_T.*log(con.Nao./con.Nai);
EK=con.E_T.*log(con.Ko./con.Ki);
EKs=con.E_T.*log((con.Ko+0.12.*con.Nao)./(con.Ki+0.12.*con.Nai));
ECa=(con.E_T./2).*log(con.Cao./Casub);

%4-aminopyridine-sensitive currents, I_toand I_sus
It0=con.g_to.*(Vm-EK).*q.*r;
Isus=con.g_sus.*(Vm-EK).*r;

%Ca2+ background current, IbCa
Ib_Ca=con.g_bCa.*(Vm-ECa);

%K+ background current, IbK
Ib_K=con.g_bK.*(Vm-EK);

%Na+ background current, IbNa
Ib_Na=con.g_bNa.*(Vm-ENa);

%L-type channel current, ICaL
ICaL12=con.g_CaL12.*(Vm-con.E_CaL).*dL12.*fL12.*fCa;
ICaL13=con.g_CaL13.*(Vm-con.E_CaL).*dL13.*fL13.*fCa;
ICaL=ICaL12+ICaL13;

%T-type Ca2+ current, ICaT
ICaT=con.g_CaT.*(Vm-con.E_CaT).*dT.*fT;

%Hyperpolarization-activated, funny current, If

IfNa=0.3833*con.g_If.*(Vm-ENa).*y;
IfK=0.6167*con.g_If.*(Vm-EK).*y;
If=IfNa+IfK;

%K1 current, IK1
xk1inf = 1./(1+exp(0.070727.*(Vm-EK)));
IK1=con.g_K1.*xk1inf.*(con.Ko./(con.Ko+0.228880)).*(Vm-EK);

%Kr current, IKr
IKr=con.g_Kr.*(Vm-EK).*pa.*pi;

%Ks current, IKs
IKs=con.g_Ks.*(Vm-EKs).*(xs.^2);

%Sodium current, INa
FNa=((9.52e-02).*exp((-6.3e-2).*(Vm+34.4))./(1+1.66.*exp(-0.225.*(Vm+63.7))))+8.69e-2;
hs11=(1-FNa).*h11+FNa.*j11;
hs15=(1-FNa).*h15+FNa.*j15;
INa11=con.g_Na11.*(m11.^3).*hs11.*Vm.*con.Nao.*con.F./(con.E_T.*1000).*(exp((Vm-ENa)./con.E_T)-1)./(exp(Vm./con.E_T)-1);
INa15=con.g_Na15.*(m15.^3).*hs15.*Vm.*con.Nao.*con.F./(con.E_T.*1000).*(exp((Vm-con.ENa15)./con.E_T)-1)./(exp(Vm./con.E_T)-1);
INa=INa11+INa15;

%Na+ - K+ pump current, INaK
INaK=con.I_NaKmax.*(((con.Ko).^1.2)./((con.K_mK).^1.2+(con.Ko).^1.2)).*((con.Nai.^1.3)./((con.K_mNa).^1.3+(con.Nai).^1.3))./(1.0+exp(-(Vm-ENa+120)/30));

%Na+-Ca2+ exchanger current, INCX

d0=1+(con.Cao./con.K_co).*(1+exp(con.Q_co.*Vm./con.E_T))+(con.Nao./con.K_1no).*(1+(con.Nao./con.K_2no).*(1+con.Nao./con.K_3no)); 
k43=con.Nai./(con.K_3ni+con.Nai);
k41=exp(-con.Q_n.*Vm./(2*con.E_T));
k34=con.Nao./(con.K_3no+con.Nao);
k21=(con.Cao./con.K_co).*exp(con.Q_co.*Vm./con.E_T)./d0;
k23=(con.Nao./con.K_1no).*(con.Nao./con.K_2no).*(1+con.Nao./con.K_3no).*exp(-con.Q_n*Vm./(2*con.E_T))./d0;
k32=exp(con.Q_n.*Vm./(2.*con.E_T));
x1=k34.*k41.*(k23+k21)+k21.*k32.*(k43+k41);
di=1+(Casub./con.K_ci).*(1+exp(-con.Q_ci.*Vm./con.E_T)+con.Nai./con.K_cni)+(con.Nai./con.K_1ni).*(1+(con.Nai./con.K_2ni).*(1+con.Nai./con.K_3ni));
k12= (Casub./con.K_ci).*exp(-con.Q_ci.*Vm./con.E_T)./di;
k14=(con.Nai./con.K_1ni).*(con.Nai./con.K_2ni).*(1+con.Nai./con.K_3ni).*exp(con.Q_n.*Vm./(2*con.E_T))./di;
x2=k43.*k32.*(k14+k12)+k41.*k12.*(k34+k32);
x3=k43.*k14.*(k23+k21)+k12.*k23.*(k43+k41);
x4=k34.*k23.*(k14+k12)+k21.*k14.*(k34+k32);

INaCa=con.K_NaCa.*(k21.*x2-k12.*x1)./(x1+x2+x3+x4);

%Sustained inward current, Ist
I_st=con.g_st.*(Vm-con.E_st).*qa.*qi;

%Acetylcholine-activated K+ current, IKACh
I_KACh=con.C*con.g_KACh_max.*(Vm-EK).*w;

%Ca2+ fluxes in the SR
jSRCarel=con.k_s.*OO.*(CajSR-Casub);

%Intracellular Ca2+ flux
fPLB=2.9102*PLB.^9.5517./(0.2763^9.5517+PLB.^9.5517)+0.4998;
j_up=(con.P_up.*fPLB)./(1+con.K_up./Cai);

j_Cadiff=(Casub-Cai)./con.tho_difCa; %Ca2+ diffusion flux from submembrane space to myoplasm, jCa,dif
j_tr=(CanSR-CajSR)./con.tho_tr; %Ca2+ flux between network and junctional SR compartments, jtr

ATP=con.ATP_max.*((con.kATP.*(cAMP*100./con.cAMPb).^con.n_ATP)./(con.k_ATP05+(cAMP*100./con.cAMPb).^con.n_ATP)-con.K_ATPmin)/100;
PKA=update_PKA(cAMP);

I=[It0 Isus Ib_Ca Ib_K Ib_Na ICaL ICaT If IK1 IKr IKs INa INaK INaCa I_st I_KACh jSRCarel j_up j_Cadiff j_tr ATP PKA];
