function [ PowerTerms, steerdirectionterms ] = CSAProdPowTermsSeq( L,Me,Ne,M,N,beta,u0 )
% PowerTerms = CSAProdPowTermsSeq( M,N,Me,Ne,L,beta )
% computes terms in sequence alpha^{|aM-bN|+|dM-eN|} and
% alpha^[M|a-d|+N|b-e|]
%#codegen
% L = M*N*beta;
% Me = L/M;
% Ne = L/N;
    PowerTerms = zeros(1,2*L^2*beta^2);
    steerdirectionterms = zeros(1,2*L^2*beta^2);
    ind = [1,2];
    for a = 0:Me-1
        for b = 0:Ne-1
            for d = 0:Me-1
                for e = 0:Ne-1;
                    PowerTerms(ind) = [abs(a*M-b*N)+abs(d*M-e*N),...
                                       M*abs(a-d)+N*abs(b-e)];
                    steerdirectionterms(ind) = [exp(-1i*pi*(a*M-b*N)*u0)*exp(1i*pi*(d*M-e*N)*u0), ...
                                                exp(-1i*pi*(a*M-d*M)*u0)*exp(1i*pi*(b*N-e*N)*u0)];
                    ind = ind+2;
                end
            end
        end
    end
    
    
end

