function K = subspace_qkd(v,k,d,m)

	p = v+(1-v)*k/d;
	K = p*mub_qkd(v/p,k,m);
	anal = p*Kiso(v/p,k);
	anal_min = p*Kiso_min(v/p,k);

end

function K = Kiso(v,d)

	vt = v + (1-v)/d^2;
	K = log2(d) - (1-vt)*log2(d^2-1) - binary_entropy(vt);
		
end

function K = Kiso_min(v,d)

	W = v+(1-v)/d;
	P = 1 - W + 1/d + (2*(-1 + W))/d^2 + ( 2*sqrt((-1 + d)*(1 - W)*(-1 + W + d*W)))/d^2;
	K = -log2(P) - hab(v,d);
	
end


