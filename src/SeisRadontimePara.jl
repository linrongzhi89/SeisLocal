"""
    LOCAL_RADON_FORWARD_ADJOING: Operators for t-x, t-q Radon transforms

[Out] = radon_forward_adjoint(In, Par, type_of_transform);

IN    Radon panel if itype =  1 (Forward transform)  with In(nq,nt)
        CMP gather  if itype = -1 (Adjoint transform)  with In(nh,nt)

OUT   CMP gather  if itype =  1 (Forward transform)  with Out(nh,nt)
        Radon panel if itype = -1 (Adjoint transform)  with Out(nq,nt)

Par.offset  :  vector containing the nh offsets
Par.offset  :  vector containing the nq curvatures
Par.dt      :  sampling interval
Par.wavelet :  zero phase wavelet or band-liminting operator with
                 bw similar to that of the data
Par.L       :  window length for local Radon transform (Total length is 2L)

itype =  1  :  forward(from Radon domian to time domian)
itype = -1  :  adjoint(from time domian to Radon domian)

Copyright (C) 2019, Signal Analysis and Imaging Group
For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
Author: M.D.Sacchi modified, by Rongzhi (Matlab to Julia)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details: http://www.gnu.org/licenses/
"""
function SeisRadontimePara(In::Array{Tin,2}, adj::Bool;
                        nt=100, dt=0.004, R=Ricker(),
                        h=collect(0.0:20.0:1000.0),
                        q=collect(-0.05:0.01:2.2))where{Tin<:Real}

    nq = length(q);
    nh = length(h);
    h2 = (h/maximum(abs.(h))).^2;

    if adj
        m = In;
        (nt,nq) = size(m);
        d = zeros(nt,nh);
         R = reshape(R,length(R),1)
        mcon = conv2(m,R);
        (nt2,nq2) = size(mcon);
        m = mcon[Int(ceil(nt2/2)-ceil(nt/2)+1):Int(ceil(nt2/2)+ceil(nt/2)),:];
    else
        d = In;
        (nt,nh) = size(d);
        m = zeros(nt,nq);
    end


    for it = 1:nt
        for iq = 1:nq
            for ih = 1:nh
                tau = (it-1)*dt - q[iq]*h2[ih];
                itau = tau/dt+1;
                itau1 = Int(floor(itau));
                itau2 = itau1 + 1;
                a= itau-itau1;
                if itau2<nt && itau1>1
                    if adj
                        d[it,ih] = d[it,ih] + (1-a) * m[itau1,iq] + a*m[itau2,iq];
                    else
                        m[itau1,iq] = m[itau1,iq] + (1-a) * d[it,ih];
                        m[itau2,iq] = m[itau2,iq] + (a) * d[it,ih];
                    end
                end

            end
        end
    end

    if adj
        Out = d;
    else
        R = reshape(R,length(R),1)
        mcon = conv2(m,R);
        (nt2,nq2) = size(mcon);
        m = mcon[Int(ceil(nt2/2)-ceil(nt/2)+1):Int(ceil(nt2/2)+ceil(nt/2)),:];
        Out = m;
    end
    return Out
end


































# function Ricker(; dt=0.002, f0=20.0)
#
#     nw = 2.0/(f0*dt)
#     nc = floor(Int, nw/2)
#     t = dt*collect(-nc:1:nc)
#     b = (pi*f0*t).^2
#     w = (1 .- 2 .*b).*exp.(-b)
# end
