function confidence_region_metropolis(d,mu,sigma,zeta,sigmastep,base)

	npars = length(mu);
	alpha = 0.05;
	zetaguess = sqrt(2*gammaincinv(1-alpha,npars/2));
	if zeta == 0
		zeta = zetaguess
	end

	sigmainv = inv(sigma);

function sample = metropolis(current)

	proposal = current + sigmastep*randn(npars,1);
	if min(proposal) < 0
		sample = current;
	elseif mineig(d,proposal,base) < 0
		sample = current;
	else
		density_proposal = exp(-(proposal-mu)'*sigmainv*(proposal-mu)/2);
		density_current = exp(-(current-mu)'*sigmainv*(current-mu)/2);
		alpha = density_proposal/density_current;
		u = rand();
		if u <= alpha
			sample = proposal;
			nacceptance = nacceptance + 1;
		else
			sample = current;
		end
	end
end

	rho0 = initial_state(d,mu,sigma,base);
	v = 0.999;
	rho0 = v*rho0 + (1-v)*eye(d^2)/d^2;
	current = real(base*rho0(:));


	ntests = 1000;
	nburn = 1000;
	nstep = 5;
	nacceptance = 0;

	for i=1:nburn
		current = metropolis(current);
	end
	%muaverage = zeros(npars,1);
	ninterval = 0;
	for i=1:ntests
		current = metropolis(current);
	%	muaverage = muaverage + current;
		if (current-mu)'*sigmainv*(current-mu) <= zeta^2
			ninterval = ninterval + 1;
		end
		for i=1:nstep
			current = metropolis(current);
		end
	end 
	
	p_acceptance = nacceptance/(nburn+ntests*(nstep+1))
	confidence = ninterval/ntests
	%muaverage = muaverage/ntests

end


function rho = initial_state(d,mu,sigma,base)

	yalmip('clear');

	rho = sdpvar(d^2,d^2,'hermitian','complex');
	delta = sdpvar;

	v = (base*rho(:) - mu);

	constraints = [rho >= 0, trace(rho) == 1, [delta v'; v sigma] >= 0];

	ops = sdpsettings(sdpsettings,'verbose',0,'solver','sedumi','cachesolvers',1);
	optimize(constraints,delta,ops);
	rho = value(rho);

end 

function lambda = mineig(d,mu,base)

	yalmip('clear');

	rho = sdpvar(d^2,d^2,'hermitian','complex');
	lambda = sdpvar;

	constraints = [rho - lambda*eye(d^2) >= 0, trace(rho) == 1, base*rho(:) == mu];

	ops = sdpsettings(sdpsettings,'verbose',0,'solver','mosek','cachesolvers',1);
	optimize(constraints,-lambda,ops);
	lambda = value(lambda);

end 
