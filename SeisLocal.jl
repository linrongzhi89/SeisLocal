"""
Function to patch 2D seismic data.
Syntax:
%        [OUT,Minval,Maxval]=SeisPatch2D(IN,dt,it_WL,it_WO,dx,ix_WL,ix_WO,xmin,xmax)
% INPUTS:
%        inp   - Full seismic section
%        dt    - sampling in time
%        it_WL - window length in time
%        it_WO - window overlap in time (give in %)
%        dx    - sampling in space
%        ix_WL - window length in space
%        ix_WO - window overlap in space (give in %)
%        xmin  - minimum number of trace
%        xmax  - maximum number of trace
%
% OUTPUTS:
%        OUT - Patches
%        Minval - vector with minimum time/space in patches
%        Maxval - vector with maximum time/space in patches
"""
function SeisLocal(inp; dt=0.004,it_WL=20,it_WO=20,dx=1,ix_WL=20,ix_WO=20,xmin=1,xmax=100)

    nt = size(inp,1);
    nx = xmax-xmin+1;

    # Window overlaps in samples
    it_WO  = round(Int64,(it_WO*it_WL)/100);
  	ix_WO  = round(Int64,(ix_WO*ix_WL)/100);

    # Safe-gard for window bigger than data
    if it_WL > nt
        it_WL = nt;
    end

    if ix_WL > nx
        ix_WL = nx;
    end

    tmax = dt*(nt-1);
    #  number of windows in t and x
    it_NW = floor(Int64,nt/(it_WL-it_WO));
    ix_NW = floor(Int64,nx/(ix_WL-ix_WO));

    # % Safe-guard to account for ovelarpping samples
    # % This assumes data starts from 0 time

    if (dt*(it_NW-1)*(it_WL-it_WO)+dt*(it_WL)) < tmax
        it_NW = it_NW+1;
    end

    if (xmin+(ix_NW-1)*(ix_WL-ix_WO)+ix_WL) < xmax
        ix_NW = ix_NW+1;
    end

    # % Number of patches

    NW=it_NW*ix_NW;
    println("Number of patches= ", NW);
    np = it_NW*ix_NW;
    Minval = zeros(np,2);
    Maxval = zeros(np,2);
    # % Split section into n multipatches
    npatch=0;

    for it_W=1:it_NW
        # % Get starting and ending points of patches
        mint=dt*(it_W-1)*(it_WL-it_WO);
        maxt=mint+dt*(it_WL-1);

        if maxt>=tmax
            maxt=tmax;
        end

        for ix_W=1:ix_NW
            minx1=xmin+(ix_W-1)*(ix_WL-ix_WO);
            maxx1=minx1+(ix_WL-1);

            if  maxx1 >= xmax
                maxx1 =  xmax;
            end
            npatch=npatch+1;
            Minval[npatch,:]=[mint minx1];
            Maxval[npatch,:]=[maxt maxx1];
        end
    end
    # % Actual patch extraction
    OUT = Array{Float64,3}(zeros(it_WL,ix_WL,npatch));
    for ipatch=1:npatch
        temp = SeisLocalPatch(inp,Minval[ipatch,:],Maxval[ipatch,:];dt=0.004);
        OUT[1:size(temp,1),1:size(temp,2),ipatch] = temp;
    end

    OUT = Array{Float64,3}(OUT);
    return OUT, Minval, Maxval
end
