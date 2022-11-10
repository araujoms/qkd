function K = mehul_experiment_restricted()

	load('mehul_data','data')

	d=3;
	counts = data{d};
	%counts = 10*data{d}; %simulated data

	dirty_counts = counts;
	dirty_counts(dirty_counts == 0) = 1;

	p = zeros(d,d,d+1);
	dirty_p = zeros(d,d,d+1);
	ncounts = zeros(d+1,1);
	ndirty_counts = zeros(d+1,1);
	for i=1:d+1
		ncounts(i) = sum(sum(counts(:,:,i)));
		p(:,:,i) = counts(:,:,i)/ncounts(i);
		ndirty_counts(i) = sum(sum(dirty_counts(:,:,i)));
		dirty_p(:,:,i) = dirty_counts(:,:,i)/ndirty_counts(i);
	end

	f_exp = zeros((d+1)*(d-1),1);
	dirtyf_exp = zeros((d+1)*(d-1),1);
	base = restricted_mub_basis(d);

	counter = 0;
	for k=1:d+1
		for l=0:d-2
			counter = counter + 1;
			for i=0:d-1
				f_exp(counter) = f_exp(counter) + p(i+1,mod(i+l,d)+1,k);
				dirtyf_exp(counter) = dirtyf_exp(counter) + dirty_p(i+1,mod(i+l,d)+1,k);
			end
		end
	end

	sigma = zeros((d+1)*(d-1),(d+1)*(d-1));
	for k=1:d+1
		s = (k-1)*(d-1)+1;
		e = k*(d-1);
		sf_exp = f_exp(s:e);
		sigma(s:e,s:e) = (-sf_exp*sf_exp' + diag(sf_exp) )/ndirty_counts(k);
	end

	alpha = 0.05;
	npars = length(f_exp);
	zetaguess = sqrt(2*gammaincinv(1-alpha,npars/2));
	zeta = 0.9*zetaguess; %with real data
	%zeta = 1.05*zetaguess; %with simulated data
	sigmastep = 0.0004; %with real data
	%sigmastep = 0.0002; %with simulated data
	%confidence_region_metropolis(d,f_exp,sigma,zeta,sigmastep,base)

	K = hae_mub_restricted(f_exp,sigma,zeta,d,8) - conditional_entropy(p(:,:,1));

end

function objective = hae_mub_restricted(f_exp,sigma,zeta,d,m)

	[w,t] = gauss_radau(m);

	basis = restricted_mub_basis(d);

	yalmip('clear');

	rho = sdpvar(d^2,d^2,'hermitian','complex');
	z = sdpvar(d^2,d^2,d,m,'full','complex');
	zdz = sdpvar(d^2,d^2,d,m,'hermitian','complex');
	zzd = sdpvar(d^2,d^2,d,m,'hermitian','complex');

	f = real(basis*rho(:)) - f_exp;

	objective = 0;
	for j=1:m
		for i=1:d
			partAB = kron(diag(ket(i,d)),eye(d));
			partE = z(:,:,i,j) + z(:,:,i,j)' + (1-t(j))*zdz(:,:,i,j);
			objective = objective +  (w(j)/(t(j)*log(2)))*(partAB(:)'*partE(:) + t(j)*trace(zzd(:,:,i,j)));
		end
	end

	constraints = [trace(rho) == 1, [zeta^2 f'; f sigma] >= 0];
	for i=1:d
		for j=1:m
			constraints = [constraints, [rho, z(:,:,i,j); z(:,:,i,j)', zdz(:,:,i,j)] >= 0, [rho, z(:,:,i,j)'; z(:,:,i,j), zzd(:,:,i,j)] >= 0];
		end
	end

    ops = sdpsettings(sdpsettings,'verbose',0,'solver','mosek');
    sol = optimize(constraints,objective,ops);
    objective = sum(w./t)/log(2) + value(objective);
    
end

function base = restricted_mub_basis(d)

	load('paper_mubs','mubs');
	mub = mubs{d};

	base = zeros((d-1)*(d+1),d^4);

	counter = 0;
	for k=1:d+1
		for l=0:d-2
			counter = counter + 1;
			for i=0:d-1
				temp = kron(ketbra(mub(:,i+1,k)),transpose(ketbra(mub(:,mod(i+l,d)+1,k))));
				base(counter,:) = base(counter,:) + temp(:).';
			end
		end
	end

end
