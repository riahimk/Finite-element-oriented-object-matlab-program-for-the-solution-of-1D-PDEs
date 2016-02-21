classdef fem_1D < handle
    properties
        nQpt;
        Qrw;
        Qrp;
        matrix;
    end
    methods
        %[ reference_w, reference_q ] = quadrature_set ( quad_num );
        function P1=fem_1D()
            P1.nQpt=2;
            set_Gauss_Legendre_quad(P1.nQpt);
            function [q_elm_posi]=maps(K_ref,ref_Q_posi,elm_posi)
                q_elm_posi = ( ( max(K_ref(:)) - ref_Q_posi(1,1) ) * elm_posi(1,1)   ...
                    + (                    ref_Q_posi(2,1) - min(K_ref(:)) ) * elm_posi(1,2) ) ...
                    / ( max(K_ref(:)) - min(K_ref(:)) ); % normalization of the distance
                return;
            end
            function set_Gauss_Legendre_quad(Quad_order)
                if ( Quad_order == 1 )
                    P1.Qrp(1) = 0.0;
                    P1.Qrw(1) = 2.0;
                elseif ( Quad_order == 2 )
                    P1.Qrp(1) = - 0.577350269189625764509148780502;
                    P1.Qrp(2) =   0.577350269189625764509148780502;
                    P1.Qrw(1) = 1.0;
                    P1.Qrw(2) = 1.0;
                elseif ( Quad_order == 3 )
                    P1.Qrp(1) = - 0.774596669241483377035853079956;
                    P1.Qrp(2) =   0.0;
                    P1.Qrp(3) =   0.774596669241483377035853079956;
                    P1.Qrw(1) = 5.0 / 9.0;
                    P1.Qrw(2) = 8.0 / 9.0;
                    P1.Qrw(3) = 5.0 / 9.0;
                elseif ( Quad_order == 4 )
                    P1.Qrp(1) = - 0.861136311594052575223946488893;
                    P1.Qrp(2) = - 0.339981043584856264802665759103;
                    P1.Qrp(3) =   0.339981043584856264802665759103;
                    P1.Qrp(4) =   0.861136311594052575223946488893;
                    P1.Qrw(1) = 0.347854845137453857373063949222;
                    P1.Qrw(2) = 0.652145154862546142626936050778;
                    P1.Qrw(3) = 0.652145154862546142626936050778;
                    P1.Qrw(4) = 0.347854845137453857373063949222;
                elseif ( Quad_order == 5 )
                    P1.Qrp(1) = - 0.906179845938663992797626878299;
                    P1.Qrp(2) = - 0.538469310105683091036314420700;
                    P1.Qrp(3) =   0.0;
                    P1.Qrp(4) =   0.538469310105683091036314420700;
                    P1.Qrp(5) =   0.906179845938663992797626878299;
                    P1.Qrw(1) = 0.236926885056189087514264040720;
                    P1.Qrw(2) = 0.478628670499366468041291514836;
                    P1.Qrw(3) = 0.568888888888888888888888888889;
                    P1.Qrw(4) = 0.478628670499366468041291514836;
                    P1.Qrw(5) = 0.236926885056189087514264040720;
                elseif ( Quad_order == 6 )
                    P1.Qrp(1) = - 0.932469514203152027812301554494;
                    P1.Qrp(2) = - 0.661209386466264513661399595020;
                    P1.Qrp(3) = - 0.238619186083196908630501721681;
                    P1.Qrp(4) =   0.238619186083196908630501721681;
                    P1.Qrp(5) =   0.661209386466264513661399595020;
                    P1.Qrp(6) =   0.932469514203152027812301554494;
                    P1.Qrw(1) = 0.171324492379170345040296142173;
                    P1.Qrw(2) = 0.360761573048138607569833513838;
                    P1.Qrw(3) = 0.467913934572691047389870343990;
                    P1.Qrw(4) = 0.467913934572691047389870343990;
                    P1.Qrw(5) = 0.360761573048138607569833513838;
                    P1.Qrw(6) = 0.171324492379170345040296142173;
                else
                    fprintf ( 1, '\n' ); fprintf ( 1, '  Please request Gauss_Quad order < %d.\n', Quad_order ); error ( ' FATAL ERROR !' );
                end
            end
        end
    end
end