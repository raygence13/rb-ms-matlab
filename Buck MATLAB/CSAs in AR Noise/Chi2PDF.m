function Z = Chi2PDF(x,var)

Z = 1/(2*var) .* exp(-(x/(2*var)));
end