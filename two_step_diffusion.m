function out = two_step_diffusion(t,D1,D2,A,C1,C2,nmax,sigma,dx1,dx2,dy,rlim,varargin)
    if numel(varargin) > 0
        var = varargin{1};
        val = varargin{2};
        if var == "t0"
            t0 = val;
        else
            error("Invalid name/value pair. Only t0 is supported.")
        end
    else
        t0 = 0;
    end
    out = diffusion_moving_beam(t,D1,A,C1,nmax,sigma,dx1,dy,'rlim',rlim,"t0",t0)...
        + diffusion_moving_beam(t,D2,A,C2,nmax,sigma,dx2,dy,'rlim',rlim,"t0",t0);
end