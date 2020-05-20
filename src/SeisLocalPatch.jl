"""
% Function to extract patches from seismic.
  % Syntax:
  %        [OUT]=SeisWindowPatch(IN,Minval,Maxval;dt)
  % INPUTS:
  %        IN - Seismic section
  %        Minval - vector with minimum time/space in patches
  %        Maxval - vector with maximum time/space in patches
  %
  % OUTPUTS:
  %        OUT - Patch
  %
  % Note: Adapted from SeismicJulia package, modified from Breno.
"""
function SeisLocalPatch(IN,Minvalue,Maxvalue;dt=0.004)

    (nt,nx)=size(IN);


     mint = Minvalue[1];
    ixmin = Int64(Minvalue[2]);


     maxt = Maxvalue[1];
    ixmax = Int64(Maxvalue[2]);


   itmin=round(Int64,mint/dt)+1;
   itmax=round(Int64,maxt/dt)+1;


   itmin = itmin < 1 ? 1 : itmin;

   # if itmax>nt
   #    itmax=nt;
   # end
   itmax = itmax > nt ? nt : itmax;

   # if ixmin<1
   #    ixmin=1;
   # end
    ixmin = ixmin < 1 ? 1 : ixmin;

   # if ixmax>nx
   #    ixmax=nx;
   # end
   ixmax = ixmax > nx ? nx : ixmax;

   OUT=zeros(floor(Int64,(itmax-itmin+1)),floor(Int64,(ixmax-ixmin+1)));

   # % Extract samples

   for ix = ixmin:ixmax
       for it = itmin:itmax
           OUT[it-itmin+1,ix-ixmin+1]=IN[it,ix];
       end
   end
   return OUT
end
