function result = is_trans_nlin_sdsos(rec1, rec2, dyn, tot_deg)
  % is_trans_nlin(rec1, rec2, dyn): Return false if a certificate that 
  % guarantees the non-existence of a flow from rec1 to rec2 under 
  % nonlinear dynamics dyn is found, true otherwise
  % 
  % Inputs 
  %   - rec1, rec2: sets
  %   - dyn = {fx, x, (d, drec)}: dynamics \dot x = fx(x,d), d \in drec
  %   - tot_deg: relaxation order for moment hierarchy

  global ops;
  result = true;

  irec = intersect(rec1,rec2);
  % overlapping
  if irec.isFullDim
    warning('isTransLinRec: was called with overlapping Recs - returning true');
    result = true;
    return;
  end

  % are they adjacent?
  if length(irec.getFlatDims) ~= 1
    result = false;
    return;
  end

  % find "flat" dimension
  flatdim = irec.getFlatDims;
  h1 = zeros(1, length(dyn{1}));
  if rec1.xmax(flatdim)==rec2.xmin(flatdim)
    h1(flatdim)=1;
  else
    h1(flatdim)=-1;
  end 
  % h1 normal vector

  n_x = length(dyn{1});

  obj = h1*dyn{1};
    
  % Convert objective to coef-mono form
  [obj_coef, mons] = coefficients(obj);
  obj_mons = zeros(n_x, length(obj_coef));
  for i=1:length(obj_coef)
    obj_mons(:,i) = degree(mons(i), dyn{2});
  end

  [r, res] = mosekopt('symbcon echo(0)'); 
  prob = []; 

  % Basic constraint:
  %  1. gamma -si * gxi - K = obj
  % 
  % Variable vector [gamma  K  s1 ... si]

  sigma_deg = tot_deg - 2;   % degree of all g's is 2

  Ag = PolyLinTrans.eye(1, n_x, 0, tot_deg).as_vector_trans;
  Ak = PolyLinTrans.eye(n_x, n_x, tot_deg, tot_deg).as_matrix_trans;

  n_g = 1;
  n_k = size(Ak, 2);

  % Get Asi's
  Asi = {};
  nsi = [];

  for i=1:n_x
    xm = irec.xmin(i);
    xp = irec.xmax(i);
    coef = [-xm*xp xm+xp -1];
    mons = zeros(n_x, 3);
    mons(i,:) = [0 1 2];
    Asi{end+1} = PolyLinTrans.mul_pol(n_x, sigma_deg, mons, coef).as_matrix_trans;
    nsi(end+1) = size(Asi{end},2);
  end

  % [Ag   -Ak   -As1 -As2] =   [obj];
  As = [Asi{:}];
  prob.a = [Ag  -Ak  -As];
  prob.blc = zeros(size(Ag, 1), 1);
  prob.buc = zeros(size(Ag, 1), 1);

  for i=1:length(obj_coef)
    prob.blc(mono_rank_grlex(n_x, obj_mons(:,i)), 1) = obj_coef(i);
    prob.buc(mono_rank_grlex(n_x, obj_mons(:,i)), 1) = obj_coef(i);
  end

  % Make K sdsos
  prob = add_sdsos(prob, 1+n_g, n_g + n_k);

  % Make the si's sdsos
  for i=1:n_x
    s1 = 1 + n_g + n_k + sum(nsi(1:i-1));
    s2 = n_g + n_k + sum(nsi(1:i));
    prob = add_sdsos(prob, s1, s2);
  end

  prob.c = zeros(1, size(prob.a, 2));
  prob.c(1) = 1;

  prob.blx = -inf*ones(1, size(prob.a,2));
  prob.bux = +inf*ones(1, size(prob.a,2));

  [r, res] = mosekopt('minimize echo(0)', prob); 

  gamma_opt = res.sol.itr.xx(1);

  if strcmp(res.sol.itr.solsta, 'OPTIMAL') && gamma_opt < 0
    result = false;
  end
end

function prob = add_sdsos(prob, s1, s2)
  % Add an sdsos constraint to 'prob' for the variables x(s1:s2)

  if size(prob.a,1) ~= length(prob.buc)
    error('dimension mismatch')
  end

  [numcon, numvar] = size(prob.a);
  
  L = 1+s2-s1;  % number of variables

  if s1+L > numvar
    error('not enough variables')
  end

  n = (sqrt(1+8*L) - 1)/2;
  if ceil(n) ~= n
    error(strcat('L=', num2str(L), ' does not represent a symmetric matrix'))
  end

  D = sdd_mat(L);

  [d1 d2] = size(D);

  % A = [ -Aold-   0 
  %      O  -I  0  D]
  % with I at pos s1:s1+L

  prob.a = [prob.a            sparse(numcon, d2); 
            sparse(d1, s1-1)  -speye(d1)  sparse(d1, numvar-d1-s1+1)  D];
  prob.buc = [prob.buc; zeros(d1, 1)];
  prob.blc = [prob.blc; zeros(d1, 1)];

  num_cones = d2/3;

  % 1 represents rotated cone 2 x1 x2 >= x3 (var res.symbcon.MSK_CT_RQUAD)
  try
    prob.cones.type = [prob.cones.type ones(1, num_cones)];
    prob.cones.subptr = [prob.cones.subptr length(prob.cones.sub)+(1:3:3*num_cones)];
    prob.cones.sub = [prob.cones.sub numvar+(1:3*num_cones)];
  catch
    prob.cones.type = ones(1, num_cones);
    prob.cones.sub = numvar+(1:3*num_cones);
    prob.cones.subptr = 1:3:3*num_cones;
  end
end


