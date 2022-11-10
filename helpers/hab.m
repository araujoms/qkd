function h = hab(v,d)
	if v==1
		h = 0;
	else
		h = binary_entropy(v+(1-v)/d) + (1-v-(1-v)/d)*log2(d-1);
	end
end
