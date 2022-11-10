function K = overlap_qkd(v,d,m)

	K = hae(v,d,m) - hab(v,d);

end

function objective = hae(v,d,m)

	[w,t] = gauss_radau(m);
	f_iso = [v/d+(1-v)/d^2;((1-v)/d^2)*ones(d,1)].*ones(d+1,5*(d-1));
	f_iso = f_iso(:);

	basis = overlap_basis(d);

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


function B = overlap_basis(d)

	localB = zeros(d,d,5);

	for i=1:d
		localB(:,i,1) = ket(i,d);
	end

	for i=1:d/2
		localB(:,2*i-1,2) = (ket(2*i-1,d) + ket(2*i,d))/sqrt(2);
		localB(:,2*i,2) = (ket(2*i-1,d) - ket(2*i,d))/sqrt(2);
	end

	for i=1:d/2
		localB(:,2*i-1,3) = (ket(2*i-1,d) + 1i*ket(2*i,d))/sqrt(2);
		localB(:,2*i,3) = (ket(2*i-1,d) - 1i*ket(2*i,d))/sqrt(2);
	end

	localB(:,1,4) = ket(1,d);
	localB(:,d,4) = ket(d,d);

	for i=1:d/2-1
		localB(:,2*i,4) = (ket(2*i,d) + ket(2*i+1,d))/sqrt(2);
		localB(:,2*i+1,4) = (ket(2*i,d) - ket(2*i+1,d))/sqrt(2);
	end

	localB(:,1,5) = ket(1,d);
	localB(:,d,5) = ket(d,d);

	for i=1:d/2-1
		localB(:,2*i,5) = (ket(2*i,d) + 1i*ket(2*i+1,d))/sqrt(2);
		localB(:,2*i+1,5) = (ket(2*i,d) - 1i*ket(2*i+1,d))/sqrt(2);
	end

	B = zeros(5*(d^2-1),d^4);

	counter = 0;
	for k=1:5
		for i=0:d-1
			for j=0:d-1
				if i ~= d-1 || j ~= d-1
					counter = counter + 1;
					temp = kron(ketbra(localB(:,i+1,k)),transpose(ketbra(localB(:,j+1,k))));
					B(counter,:) = B(counter,:) + temp(:).';
				end
			end
		end
	end

end
