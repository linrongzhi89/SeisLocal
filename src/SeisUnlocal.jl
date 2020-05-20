"""
%  (Patches,nt,dt,it_WL,it_WO,Mint,min_ix,max_ix,dx,ix_WL,ix_WO,Minx,Maxx)
   % Function to unpatch 2D seismic data.
   % Syntax:
   %        [OUT]= SeisUnlocal
   %     (Patches,Minval,Maxval,nt;dt,it_WL,it_WO,dx,ix_WL,ix_WO,xmin,xmax)
   % INPUTS:
   %        Patches - Patches of full seismic section
   %        Minval 	- vector with minimum time/space in patches
   %        Maxval 	- vector with maximum time/space in patches
   %        nt 		- number of time samples in full section
   %        dt 		- sampling in time
   %        it_WL 	- window length in time
   %        it_WO 	- window overlap in time (give in %)
   %        dx 		- sampling in space
   %        ix_WL	- window length in space
   %        it_WO 	- window overlap in space (give in %)
   %        xmin 	- minimum number of trace
   %        xmax 	- maximum number of trace
   % OUTPUTS:
   %        OUT - Unpatched (full) seismic section
"""
function SeisUnlocal(Patches,Minval,Maxval,nt;dt=0.004,it_WL=20,it_WO=20,dx=1,ix_WL=20,ix_WO=20,xmin=1,xmax=100)


      # % Patches dimensions
		nx = Int(xmax-xmin+1);


      # % Window overlaps in samples
		it_WO=round(Int64,(it_WO*it_WL)/100);
		ix_WO=round(Int64,(ix_WO*ix_WL)/100);


		# % Safe-guard for window bigger than data
		# if it_WL > nt
		# 	it_WL=nt;
		# end
        it_WL = it_WL > nt ? nt : it_WL;

		# if ix_WL > nx
		# 	ix_WL=nx;
		# end
        ix_WL = ix_WL > nx ? nx : ix_WL;

	   # % Allocate output
	   # % TO DO
	   # % Input can be a pure or full quaternion. Can generalize
	   # % this section of the code using str2func. if condition
	   # % will be necessary to differ pure and full quaternions.

      POUT=zeros(nt,nx);


      # println("Apply Tapering...")

      for ipatch=1:size(Patches,ndims(Patches))

      # % Set initial tapers
            tapti=0; taptf=0;
            tapxi=0; tapxf=0;

      # % Extract i-th patch
            d_patch = Patches[:,:,ipatch];
	        dims_patch = size(d_patch);

	        ot_patch = Minval[ipatch,1];

    	   min_it_patch = floor(Int64,ot_patch/dt)+1;
	       max_it_patch = floor(Int64,(ot_patch/dt)+dims_patch[1]);

           if max_it_patch > nt
		      max_it_patch = nt;
		   end

	       min_ix_patch = Int(Minval[ipatch,2]);
	       max_ix_patch = Int(Maxval[ipatch,2]);


           if min_it_patch > 1
               tapti = it_WO;
           end

           if max_it_patch < nt
               taptf = it_WO;
           end

           if min_ix_patch > xmin
               tapxi = ix_WO;
           end

           if max_ix_patch < xmax
               tapxf = ix_WO;
           end


        d_patch = ApplyTaper(d_patch,max_it_patch-min_it_patch+1,max_ix_patch-min_ix_patch+1,tapti,taptf,tapxi,tapxf);

        POUT[min_it_patch:max_it_patch,min_ix_patch:max_ix_patch] = POUT[min_it_patch:max_it_patch,min_ix_patch:max_ix_patch] .+ d_patch[1:length(min_it_patch:max_it_patch),1:length(min_ix_patch:max_ix_patch)];
        POUT = copy(POUT);
      # OUT=squeeze(OUT);


    end
    return POUT
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
