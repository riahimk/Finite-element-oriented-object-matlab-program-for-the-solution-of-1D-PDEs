classdef varf <  fem_1D 
    properties
        kappa;
        initial_condition;
        boundary_condition;
        Mass;
        D_xx;
        Id;
        A;
    end
    methods
        function varf=varf(Th,N)
            if(nargin<2)
                N=1;
            end
            %if(N==1 & nargin(F)==2)
            %    error('Please check the input of the class var ! \n ~ It seems you have to add the total time micro-steps ~ \n');
            %end
            P1=fem_1D();
            varf.kappa=1.0;
            Assemble_left_and_right_handside(Th);
            function Assemble_left_and_right_handside(Mh)
                varf.D_xx= zeros(numel(Mh.vrtx));
                varf.Mass = zeros(numel(Mh.vrtx));
                for elm=1:Mh.Nelm
                    elm_length         = abs(Mh.K(2,elm)-Mh.K(1,elm));
                    ref_length         = abs(Mh.K_ref(2,1)-Mh.K_ref(1,1));
                    Quad_elm_weigth    = (elm_length./ref_length) * varf.Qrw();
                    Quad_elm_posi      = maps(Mh,elm,varf.Qrp);
                    for i=1:2
                        num_elm_vrtx_i = Mh.K_vrtx_num(i,elm);
                        [vhi,dvhi]=set_basis_func(Mh.K(i,elm),Mh.K(1,elm),Mh.K(2,elm),Quad_elm_posi());
                        %src_at_Q=eval_volume_src(Quad_elm_posi()); src_at_elm = sum(src_at_Q);
                        %varf.rhs(num_elm_vrtx_i,1)=varf.rhs(num_elm_vrtx_i,1) + src_at_elm;
                        for j=1:2
                            num_elm_node_j = Mh.K_vrtx_num(j,elm);
                            [vhj,dvhj]=set_basis_func(Mh.K(j,elm),Mh.K(1,elm),Mh.K(2,elm),Quad_elm_posi());
                            vhivhj =Quad_elm_weigth.*vhj.*vhi; m=sum(vhivhj);
                            varf.Mass (num_elm_vrtx_i,num_elm_node_j) = varf.Mass (num_elm_vrtx_i,num_elm_node_j) + m;
                            dvhivhj=(varf.kappa*Quad_elm_weigth).*(dvhj.*dvhi);  s=sum(dvhivhj);
                            varf.D_xx(num_elm_vrtx_i,num_elm_node_j) = varf.D_xx(num_elm_vrtx_i,num_elm_node_j)+s;
                        end
                    end
                end
%                 if(nargin(F)==1)
%                     varf.f = varf.Mass*F(Th.vrtx);
%                     varf.f(1,1)=varf.g(1,1);
%                     varf.f(Th.nvtx,1)=varf.g(2,1);
%                     elseif(nargin(F)==2)
%                         varf.f = varf.Mass*F(:,Th.vrtx);
%                         error('');
%                 end                    
                    
                    varf.D_xx(1,:)=0.0;      varf.D_xx(1,1) = 1.0;
                    varf.D_xx(Th.nvtx,:)=0.0;varf.D_xx(Th.nvtx,Th.nvtx)=1.0;
                    
                    varf.Mass(1,:)      =0.0;      varf.Mass(1,1)            = 1.0;
                    varf.Mass(Th.nvtx,:)=0.0;      varf.Mass(Th.nvtx,Th.nvtx)=1.0;
                    varf.Id             = varf.Mass;
            end
            function [v,v_x]      = set_basis_func(Base_index_ij,x_i,x_j,x)
                if(Base_index_ij == x_j)
                    v = (x-x_i)./ (x_j-x_i);
                    v_x = 1.0  ./ (x_j-x_i);
                elseif(Base_index_ij == x_i)
                    v  = (x_j-x) ./ (x_j-x_i);
                    v_x= -1.0    ./ (x_j-x_i);
                else
                    error('Not assigned output !')
                end
                return;
            end
            function set_bc(M,g)
                if(size(M,2)==1)
                    M(1,1)=g(1,1);
                    M(Th.nvtx,1)=g(2,1);
                else
                M(1,:)=0.0;M(1,1) = 1.0;
                M(Th.nvtx,:)=0.0;M(Th.nvtx,Th.nvtx)=1.0;
                end
                return;
            end
            function [r]=integral_Q(f,q_order,quad_posi,quad_weight)
                tmp=[];
                for q=1:q_order
                    tmp(q)=quad_weight(q)*f(quad_posi);
                end
             r = sum(tmp);
            end
            function [posi_in_elm]=maps(Mh,elm,posi_in_ref)
               % posi_in_elm = ( ( Mh.K_ref(2,1) - posi_in_ref) * Mh.K(1,elm) + (posi_in_ref - Mh.K_ref(1,1))*Mh.K(2,elm) )./abs(Mh.K_ref(2,1)-Mh.K_ref(1,1));
                posi_in_elm = .5*( ( Mh.K(2,elm) - Mh.K(1,elm)) * posi_in_ref + ( Mh.K(2,elm) + Mh.K(1,elm)));
            end
            function [y]=eval_volume_src(quad_points)
                y=F(quad_points);
                return;
            end
        end
    end
end
