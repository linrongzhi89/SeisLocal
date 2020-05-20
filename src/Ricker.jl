function Ricker(; dt=0.002, f0=20.0)

    nw = 2.0/(f0*dt)
    nc = floor(Int, nw/2)
    t = dt*collect(-nc:1:nc)
    b = (pi*f0*t).^2
    w = (1 .- 2 .*b).*exp.(-b)

end
