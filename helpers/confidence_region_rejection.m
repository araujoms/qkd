function confidence_region_rejection(d,mu,sigma,zeta,base)

	npars = length(mu);
	alpha = 0.05;
	zetaguess = sqrt(2*gammaincinv(1-alpha,npars/2));
	if zeta == 0
		zeta = zetaguess
	end

	ntests = 10000;
	nvalid = 0;
	ninterval = 0;
	sigmachol = chol(sigma,'lower');
	sigmainv = inv(sigma);

	for i=1:ntests
		sample = mu + sigmachol*randn(npars,1);
		if min(sample) >= 0
			if mineig(d,sample,base) >= 0
				nvalid = nvalid + 1;
				if (sample-mu)'*sigmainv*(sample-mu) <= zeta^2
					ninterval = ninterval + 1;
				end
			end
		end
	end 

	p_valid = nvalid/ntests
	confidence = ninterval/nvalid

	end

function [lambda rho] = mineig(d,mu,base)

	yalmip('clear');
	rho = sdpvar(d^2,d^2,'hermitian','complex');
	lambda = sdpvar;

	constraints = [rho - lambda*eye(d^2) >= 0, trace(rho) == 1, real(base*rho(:)) == mu];

	ops = sdpsettings(sdpsettings,'verbose',0,'solver','mosek','cachesolvers',1);
	optimize(constraints,-lambda,ops);
	lambda = value(lambda);

end 
