function output = simulate_with_growth(params,reaction_fn,t_save)


%Load parameters
growth_total = params.total_growth_rate;
growth_tip_fraction = params.fraction_of_growth_at_tip;
lL = params.left_boundary_position;
lF = params.length_patterning_zone;
L = params.domain_length;
l0 = params.initial_digit_length;
t0 = params.time_to_initalize_pattern;
tfin = 1e4; %simulation often stops before this once digit growth reaches end of domain
frzvars = params.committed_species;
N = params.number_gridpoints;
dt = params.time_step;
numvar = length(params.diffusion_constants);
noiseIC = 0.3;
dL = L/N;
growth_uniform = (1-growth_tip_fraction)*growth_total;
growth_tip = growth_tip_fraction*growth_total;


%variables to record uniform and tip growth
DL_uniform = 0;
DL_tip = 0;

%indices to mark left hand position, total length, right hand position, and
%committed zone of digit
l_left = round(lL*N/L);
l_digit = round(l0*N/L);
l_right = 1+l_left+l_digit;
l_committed = round(lF*N/L);


%initial conditions
rng(params.random_seed);
m = zeros(N,numvar);
m = repmat(params.initial_conditions,[N 1]).*normrnd(1,noiseIC,[N 1]);
md = m;

%coefficients of diffusion
coeffs = zeros(N,numvar);
for k = 1:N
    coeffs(k,:) = -params.diffusion_constants.*pi^2.*((((k-1)/N)*N/L)^2);
end

%%setup saving variables
n_save = length(t_save);
output.m = zeros(N,numvar,1);
output.l_digit = zeros(1,1);
output.l_committed = zeros(1,1);
output.t = zeros(1,1);

%%Simulation
committed_zone = zeros(N,1);
t = 0;
i_step = 0;
i_save = 1;
if (n_save > 0)
    t_thresh = t_save(i_save);
else
    t_thresh = 1e9;
end


while t < tfin
    t = t + dt;
    i_step = i_step + 1;
    l_right = 1+l_left+l_digit;
    
    if (l_right >= (N-1))
        break
    end
    
    m_old = m;
    digit_mask = ones(N,1);
    digit_mask(1:l_left) = 0;
    digit_mask(l_right:N) = 0;
    
    %diffuse
    for ii = 1:numvar
        md(:,ii) = real(idct(dct(m(:,ii))./((1-coeffs(:,ii).*dt))));
    end
    
    %react
    m_react = dt .* cell2mat(arrayfun(@(ROWIDX) reaction_fn(md(ROWIDX,:),params.reaction_parameters)', (1:size(md,1)).', 'UniformOutput', false));
    
    %update values from reaction term
    if(params.boundary.reflection == 1)
        m = md + m_react;
        %Use ghost cells to enforce reflective boundary conditions
        for var = 1:numvar
            m(l_right:N,var) = m((l_right-1),var);
            m(1:l_left,var) = m((l_left+1),var);
        end
    else
        %Use different reaction terms outside digit
        m_outside = dt .* cell2mat(arrayfun(@(ROWIDX) params.boundary.function(md(ROWIDX,:),params.boundary.parameters)', (1:size(md,1)).', 'UniformOutput', false));
        m = md + digit_mask .* m_react + (1-digit_mask).*m_outside;
    end
    
    %freeze committed cells
    if (t > t0)
        for ii = 1:length(frzvars)
            m(:,frzvars(ii)) = committed_zone.*m_old(:,frzvars(ii)) + (1-committed_zone).*m(:,frzvars(ii));
        end
    end 
    
    %Model 2-independent growth processes:
    if (t > t0)
        committed_zone = zeros(N,1);
        committed_zone((l_left+1):(1+l_left+l_digit - l_committed)) = 1;
        %note, here we define l_F relative to l_R
        
        %1: growth at tip - simply extend domain
        DL_tip = DL_tip + growth_tip*dt/dL;
        if(DL_tip > 1)
            l_digit = l_digit + 1;
            DL_tip = 0;
        end
        
        %2: growth throught digit - extend AND stretch domain
        DL_uniform = DL_uniform + growth_uniform*dt/dL;
        if(DL_uniform > 1)
            l_right = 1+l_left+l_digit;
            m1 = m(1:l_left,:);
            m2 = m((l_left+1):(l_left+l_digit),:);
            m2r = imresize(m2,[(size(m2,1)+1) (size(m2,2))]);
            m3 = m(l_right:(N-1),:);
            m = [m1; m2r; m3];
            l_digit = l_digit + 1;
            DL_uniform = 0;
        end
    end
    
    if (t > t_thresh)
        output.m(:,:,i_save) = m;
        output.l_digit(i_save) = l_digit;
        output.l_committed(i_save) = l_committed;
        output.t(i_save) = t;
        i_save = i_save + 1;
        if (i_save > n_save)
            t_thresh = 1e9;
        else
            t_thresh = t_save(i_save);
        end
    end
    
    
    

end

output.m(:,:,i_save) = m;
output.l_digit(i_save) = l_digit;
output.l_committed(i_save) = l_committed;
output.t(i_save) = t;

output.params = params;



end
