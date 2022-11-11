function K = mub_qkd(v,d,m)
%computes the key rate for the MUB protocol with the isotropic state
%of dimension d, with visibility v, and gauss-radau quadrature level m


	K = hae(v,d,m) - hab(v,d);
	vt = v + (1-v)/d^2;
	Kiso = log2(d) - (1-vt)*log2(d^2-1) - binary_entropy(vt);

end

function objective = hae(v,d,m)

	[w,t] = gauss_radau(m);
	f_iso = [v/d+(1-v)/d^2;((1-v)/d^2)*ones(d,1)].*ones(d+1,(d+1)*(d-1));
	f_iso = f_iso(:);

	basis = mub_basis(d);

	yalmip('clear');

	rho = sdpvar(d^2,d^2,'hermitian','complex');
	z = sdpvar(d^2,d^2,d,m,'full','complex');
	zdz = sdpvar(d^2,d^2,d,m,'hermitian','complex');
	zzd = sdpvar(d^2,d^2,d,m,'hermitian','complex');

	f = real(basis*rho(:));

	objective = 0;
	for j=1:m
		for i=1:d
			partAB = kron(diag(ket(i,d)),eye(d));
			partE = z(:,:,i,j) + z(:,:,i,j)' + (1-t(j))*zdz(:,:,i,j);
			objective = objective +  (w(j)/(t(j)*log(2)))*(partAB(:)'*partE(:) + t(j)*trace(zzd(:,:,i,j)));
		end
	end

	constraints = [trace(rho) == 1, f == f_iso];
	for i=1:d
		for j=1:m
			constraints = [constraints, [rho, z(:,:,i,j); z(:,:,i,j)', zdz(:,:,i,j)] >= 0, [rho, z(:,:,i,j)'; z(:,:,i,j), zzd(:,:,i,j)] >= 0];
		end
	end

    ops = sdpsettings(sdpsettings,'verbose',0,'solver','mosek');
    sol = optimize(constraints,objective,ops);
    objective = sum(w./t)/log(2) + value(objective);
    
end

function base = mub_basis(d)

	load('paper_mubs','mubs');
	mub = mubs{d};

	base = zeros((d^2-1)*(d+1),d^4);

	counter = 0;
	for k=1:d+1
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
