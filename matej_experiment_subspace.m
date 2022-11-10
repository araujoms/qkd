function Ksub = subspace_matej_experiment()

	load('matej_data','counts')

	p = zeros(4,4,2);
	ncounts = zeros(2,1);
	for i=1:2
		ncounts(i) = sum(sum(counts(:,:,i)));
		p(:,:,i) = counts(:,:,i)/ncounts(i);
	end

	psub = zeros(2,2);
	for base=1:2
		psub(base,:) = [sum(sum(p(1:2,1:2,base)));sum(sum(p(3:4,3:4,base)))];
	end
	average_psub = sum(psub,1)/2;

	conditional_p = zeros(2,2,2,2);
	for l=1:2
		for base=1:2
			conditional_p(:,:,base,l) = p(2*l-1:2*l,2*l-1:2*l,base)/psub(base,l);
		end
	end

	f_exp = zeros(3*2,2);
	counter = 0;
	for base=1:2
		for i=0:1
			for j=0:1
				if i ~= 1 || j ~= 1
					counter = counter + 1;
					f_exp(counter,1) = conditional_p(i+1,j+1,base,1);
					f_exp(counter,2) = conditional_p(i+1,j+1,base,2);
				end
			end
		end
	end


	ncountsub = zeros(2,2);
	for i=1:2
		ncountsub(i,1) = sum(sum(counts(1:2,1:2,i)));
		ncountsub(i,2) = sum(sum(counts(3:4,3:4,i)));
	end

	sigma = zeros(6,6,2);

	for l=1:2
		sigma(:,:,l) = blkdiag((diag(f_exp(1:3,l)) -f_exp(1:3,l)*f_exp(1:3,l)')/ncountsub(1,l),(diag(f_exp(4:6,l)) -f_exp(4:6,l)*f_exp(4:6,l)')/ncountsub(2,l));
	end

	d=2;
	base = mub_basis(d);

	npars = length(f_exp(:,1));
	alpha = 0.05;
	zeta = ones(2,1)*sqrt(2*gammaincinv(1-alpha,npars/2));
	%l=2;
	%confidence_region_rejection(d,f_exp(:,l),sigma(:,:,l),zeta(l),base)

	Ksub=0;
	for l=1:2
		Ksub = Ksub + average_psub(l)*(hae_mub(f_exp(:,l),sigma(:,:,l),zeta(l),d,8) - conditional_entropy(conditional_p(:,:,1,l)));
	end

	
end

function hae = hae_mub(f_exp,sigma,zeta,d,m)

	[w,t] = gauss_radau(m);
	basis = mub_basis(d);

	yalmip('clear');

	rho = sdpvar(d^2,d^2);
	z = sdpvar(d^2,d^2,d,m,'full');
	zdz = sdpvar(d^2,d^2,d,m);
	zzd = sdpvar(d^2,d^2,d,m);

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
    hae = sum(w./t)/log(2) + value(objective);
    
end

function base = mub_basis(d)

	load('paper_mubs','mubs');
	mub = mubs{d};

	base = zeros((d^2-1)*(d),d^4);

	counter = 0;
	for k=1:d
		for i=0:d-1
			for j=0:d-1
				if i ~= d-1 || j ~= d-1
					counter = counter + 1;
					temp = kron(ketbra(mub(:,i+1,k)),transpose(ketbra(mub(:,j+1,k))));
					base(counter,:) = temp(:).';
				end
			end
		end
	end

end
