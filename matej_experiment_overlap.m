function Koverlap = overlap_matej_experiment()

	load('matej_data','counts')

	d=4;
	p = zeros(d,d,3);
	ncounts = zeros(3,1);
	for i=1:3
		ncounts(i) = sum(sum(counts(:,:,i)));
		p(:,:,i) = counts(:,:,i)/ncounts(i);
	end

	f_exp = zeros(30,1);
	counter = 0;
	for k=1:1
		for i=0:d-1
			for j=0:d-1
				if i ~= d-1 || j ~= d-1
					counter = counter + 1;
					f_exp(counter) = p(i+1,j+1,k);
				end
			end
		end
	end

	for k=2:2
		for i=0:d-1
			for j=0:d-1
				if ~all([i j] == [0 1]) && ~all([i j] == [0 3]) && ~all([i j] == [2 1]) && ~all([i j] == [2 3])
					counter = counter + 1;
					f_exp(counter) = p(i+1,j+1,k);
				end
			end
		end
	end

	k=3;
	for i=0:d-1
		for j=0:d-1
			if (i ~= 0 && i ~= d-1) && ( j ~= 0 && j ~= d-1) && ~all([i,j] == [1 2])
				counter = counter + 1;
				f_exp(counter) = p(i+1,j+1,k);
			end
		end
	end

	f1 = f_exp(1:15);
	f2 = f_exp(16:27);
	f3 = f_exp(28:30);

	sigma1 = (diag(f1)-f1*f1')/ncounts(1);
	sigma2 = (diag(f2)-f2*f2')/ncounts(2);
	sigma3 = (diag(f3)-f3*f3')/ncounts(3);

	sigma = blkdiag(sigma1,sigma2,sigma3);

	base = overlap_basis(d);

	npars = length(f_exp);
	alpha = 0.05;
	zetaguess = sqrt(2*gammaincinv(1-alpha,npars/2));
	zeta = zetaguess;
	%confidence_region_rejection(d,f_exp,sigma,zeta,base)

	Koverlap = hae_overlap(f_exp,sigma,zeta,d,8) - conditional_entropy(p(:,:,1));

end

function hae = hae_overlap(f_exp,sigma,zeta,d,m)

	[w,t] = gauss_radau(m);

	basis = overlap_basis(d);

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

function B = overlap_basis(d)

	localB = zeros(d,d,3);

	for i=1:d
		localB(:,i,1) = ket(i,d);
	end

	for i=1:d/2
		localB(:,2*i-1,2) = (ket(2*i-1,d) + ket(2*i,d))/sqrt(2);
		localB(:,2*i,2) = (ket(2*i-1,d) - ket(2*i,d))/sqrt(2);
	end

	localB(:,1,3) = ket(1,d);
	localB(:,d,3) = ket(d,d);

	for i=1:d/2-1
		localB(:,2*i,3) = (ket(2*i,d) + ket(2*i+1,d))/sqrt(2);
		localB(:,2*i+1,3) = (ket(2*i,d) - ket(2*i+1,d))/sqrt(2);
	end

	B = zeros(30,d^4);

	counter = 0;
	for k=1:1
		for i=0:d-1
			for j=0:d-1
				if i ~= d-1 || j ~= d-1
					counter = counter + 1;
					temp = kron(ketbra(localB(:,i+1,k)),transpose(ketbra(localB(:,j+1,k))));
					B(counter,:) = temp(:).';
				end
			end
		end
	end

	for k=2:2
		for i=0:d-1
			for j=0:d-1
				if ~all([i j] == [0 1]) && ~all([i j] == [0 3]) && ~all([i j] == [2 1]) && ~all([i j] == [2 3])
					counter = counter + 1;
					temp = kron(ketbra(localB(:,i+1,k)),transpose(ketbra(localB(:,j+1,k))));
					B(counter,:) = temp(:).';
				end
			end
		end
	end

	k=3;
	for i=0:d-1
		for j=0:d-1
			if (i ~= 0 && i ~= d-1) && ( j ~= 0 && j ~= d-1) && ~all([i,j] == [1 2])
				counter = counter + 1;
				temp = kron(ketbra(localB(:,i+1,k)),transpose(ketbra(localB(:,j+1,k))));
				B(counter,:) = temp(:).';
			end
		end
	end

end
