function h = binary_entropy(x)

function l = goodlog(s)
	if s == 0
		l = 0;
	else
		l = log2(s);
	end
end

h = -x*goodlog(x) - (1-x)*goodlog(1-x);

end
