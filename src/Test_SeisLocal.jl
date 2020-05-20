
using SeisPlot
using LinearAlgebra
using LocalRadon
using Random
using DSP


nq = 100

nh = 80
dt = 0.004
h = collect(range(0.0;stop=1000.0,length=nh))      
p = collect(range(0.1,stop=0.4, length=nq))
q = p
w = Ricker(dt=0.004,f0=20)
nt = 250;

m = zeros(nt,nq);
m[80,20]  = 1;
m[120,50] =-1;
m[180,80] = 1;
# m[200,85] =-1;
w2 = reshape(w,length(w),1);
m2con = conv2(m,w2)
(nt2,nq2) = size(m2con);
m2con = m2con[Int(ceil(nt2/2)-ceil(nt/2)+1):Int(ceil(nt2/2)+ceil(nt/2)),:];

param = Dict(:nt=>nt,:dt=>dt,:R=>w,:h=>h,:q=>p);
d0 = SeisRadontimePara(m, true; param...);


sigma = 0.08 *maximum(d0);
d = d0 .+ sigma*randn(size(d0))


d = randn(size(d));

close("all")


param2 = Dict(:dt=>0.004,:it_WL=>50,:it_WO=>20,:dx=>1,:ix_WL=>20,:ix_WO=>20,:xmin=>1,:xmax=>80);


(OUT,Minval,Maxval) = SeisLocal(d; param2...)


P = SeisUnlocal(OUT,Minval,Maxval,nt; param2...)

println("error ",norm(P-d))

SeisPlotTX(OUT[:,:,11], title="Local Data", fignum=2,xlabel="Offset [m]",style="color",vmin=-1,vmax=1)

SeisPlotTX(P, title="Recover Data", fignum=3,xlabel="Offset [m]",style="color",vmin=-1,vmax=1)

SeisPlotTX(P-d, title="Difference", fignum=4,xlabel="Offset [m]",style="color",vmin=-1,vmax=1)
