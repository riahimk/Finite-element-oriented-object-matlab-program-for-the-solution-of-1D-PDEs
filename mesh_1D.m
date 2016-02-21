classdef mesh_1D < handle
    properties
        nvtx;
        vrtx;
        Nelm;
        K;
        K_vrtx_num
        h;
        K_ref;
    end
    methods
        function Th=mesh_1D(a,b,h)
            Th.h=h;
            Th.nvtx=floor(abs(b-a)/h);
            Th.Nelm=Th.nvtx-1;
            K_ref=[-1,1];
            set_vrtx(Th.nvtx,a,b);
            set_elements(Th.Nelm,Th.vrtx);
            Th.K_ref=[-1;1];
            function set_vrtx(n,a,b)
                x_max=max(a,b);x_min=min(a,b);
                Th.h = ( x_max - x_min ) / ( n - 1 );
                Th.vrtx = ( linspace ( x_min, x_max, Th.nvtx ) )';
            end
            function set_elements(Nelm,vrtx)
                Th.K(1,1:Nelm) = vrtx(1:Nelm);
                Th.K(2,1:Nelm) = vrtx(2:Nelm+1);
                Th.K_vrtx_num(1,1:Nelm)=1:numel(vrtx)-1;
                Th.K_vrtx_num(2,1:Nelm)=2:numel(vrtx);
            end
        end
    end
end