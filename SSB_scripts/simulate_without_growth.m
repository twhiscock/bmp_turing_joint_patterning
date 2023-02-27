function output = simulate_without_growth(params,reaction_fn)

%Load parameters
dt = params.time_step;
L = params.domain_length;
T = params.end_time;
seed = params.random_seed;
noiseIC = 0.3;
N = params.number_gridpoints;
n = length(params.diffusion_constants);

%Initial conditions
rng(seed);
m = repmat(params.initial_conditions,[N 1]).*normrnd(1,noiseIC,[N 1]);
md = zeros(N,n);

%coefficients of diffusion
coeffs = zeros(N,n);
for k = 1:N
    coeffs(k,:) = -params.diffusion_constants.*pi^2.*((((k-1)/N)*N/L)^2);
end

%Perform simulations
t = 0;
while t < T
    t = t + dt;
    %DIFFUSE
    for ii = 1:n
        md(:,ii) = real(idct(dct(m(:,ii))./((1-coeffs(:,ii).*dt))));
    end
    %REACT
    m = md + dt .* cell2mat(arrayfun(@(ROWIDX) reaction_fn(md(ROWIDX,:),params.reaction_parameters)', (1:size(md,1)).', 'UniformOutput', false));
end

output.m = m;
output.params  = params;

end