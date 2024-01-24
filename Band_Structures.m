#-------------------------------------------------------
# BAND STRUCTURE OF SEMICONDUCTORS 
#-------------------------------------------------------
source('mystartdefaults.m')
dispersion_relation = true; compute_dos =  false;
source('CohenBergstresser1966.m')

tic
recipunit =  1.0e+10;
ekinunit = ((hbar*recipunit)^2/(2*elm))/qel;
m=1; # Choose the semiconductor (for m = 1 to 14 )
semiconductor = material_list(m,1:4);

[qpath,tix,til]=BZpath(BZstep,qs,qe,qs_str,qe_str);

#Calculating Unit vectors of reciprocal lattice 
g = zeros(4,3);
g(1:3,1) = cross(a(:,2),a(:,3))/cell_volume;
g(1:3,2) = cross(a(:,3),a(:,1))/cell_volume;
g(1:3,3) = cross(a(:,1),a(:,2))/cell_volume;

for i=1:3;
  g(4,i) = g(1:3,i)'*g(1:3,i);
end


#Build reciprocal lattice
min_norm = sqrt(min(g(4,:)));
nstep = ceil(sqrt(sqrt(cutoff))/min_norm);
fprintf('Cutoff condition requires %d positive steps along reciprocal unit lattice vectors\n',nstep);

nodes = (2 * nstep + 1)^3;
fprintf('Generate (2* %1d + 1)^3 = %4d reciprocal lattice vectors\n',nstep,nodes);

G = zeros(5,nodes);
i = 0;
for h=-nstep:nstep;
  for k=-nstep:nstep;
    for l=-nstep:nstep;
      i+=1;
      G(1:3,i) = h*g(1:3,1) + k*g(1:3,2) + l*g(1:3,3);
      G(4,i) = G(1:3,i)'*G(1:3,i);
      G(5,i) = sqrt(G(4,i));
    endfor
  endfor
end

#Sorting the G vectors according to the norm
[G(5,:),perm] = sort(G(5,:));
G(1:4,:) = G(1:4,perm);

#Finding number of G points such that |G|^2<Gs_max
n = 1;
for i=2:length(G);
  if G(4,i)<=cutoff;
    n++;
  endif
end

ngx = ceil(n/2);
n = 2*ngx-1;   #size of hamiltonian

#Defining the pseudopotential from form factors of the semiconductor
V = zeros(1,n);
a = ls(1,m);               #Lattice spacing

fprintf('\n            G(1)            G(2)            G(3)            G(4)            V\n')

ff(:,1) = [-0.7704,-0.6953,-0.5025,-0.6771,-0.6523,-0.5103,-0.5631,-0.5117,-0.5253,-0.4457,-0.4667,-0.4491,-0.3921,-0.3112]';
%ff(m,1)' % for adjusting zero of the energy scale
for i=1:n;
  sym = 0;
  asym = 0;
  if G(4,i)<=Gs_max;
    if G(4,i)==0;
      sym = ff(m,1) * Rydberg;
    endif
    if G(4,i)==3;
      sym = ff(m,2) * Rydberg;
      asym = ff(m,5) * Rydberg;
    endif
    if G(4,i)==4;
      asym = ff(m,6) * Rydberg;
    endif
    if G(4,i)==8;
      sym = ff(m,3) * Rydberg;
    endif
    if G(4,i)==11;
      asym = ff(m,7) * Rydberg;
      sym = ff(m,4) * Rydberg;
    endif
  endif
  V(1,i) = sym * cos(2*pi*G(1:3,i)'*tau(1:3,1)) - 1i*asym * sin(2*pi*G(1:3,i)'*tau(1:3,2));
  %V(1,i) =0;
  fprintf('%15.6G %15.6G %15.6G %15.6G %15.6G\n',[G(1,i),G(2,i),G(3,i),G(4,i),V(1,i)]);
end


#Defining Hamiltonian matrix
H = zeros(n,n);
#Adding potential terms in Hamiltonian
G_diff=zeros(5,1);

for j= 1:n
  for i= 1:n
    G_diff(1:3) = G(1:3,i)-G(1:3,j);
    G_diff(5) = G_diff(1:3)'*G_diff(1:3);
    if(G_diff(5)<=Gs_max)
      for k=1:n
        if(norm(G_diff(1:3)-G(1:3,k))< tol)
          H(i,j) = V(1,k);
        end
      end
    end
  end
end

#Adding Kinetic energy terms
for k=1:length(qpath)
  for i=1:n
    H(i,i) = V(1) + ekinunit*((2*pi/a)^2) * ((qpath(1:3,k)-G(1:3,i))' * (qpath(1:3,k)-G(1:3,i)));
  endfor

  if(!ishermitian(H,tol))
    printf('\nHamiltonian matrix not hermitian : fatal error.\n');
    return;
  else
    [v,ev] = eig(H);
    E = real(diag(ev));
    [E,perm] = sort(E);
    v = v(:,perm);

    for g=1:nband
      bandstructure(g,k) = E(g);
    endfor
  endif
end


% Plot the band structure 
graph_lines = plot(qpath(5, :), bandstructure, '-', 'LineWidth', 0.7, 'Color', 'k');
hold on;
graph_markers = plot(qpath(5, 1:5:end), bandstructure(:, 1:5:end), 'o', 'MarkerSize', 2.5, 'Color', 'k', 'MarkerFaceColor','w');

hold off;

xlim([0, max(qpath(5, :))]);
ylim([-4, 8]);
set(gca, 'ytick', [-16:1:24], 'FontSize', 8, 'FontWeight', 'bold');
set(gca, 'xtick', tix);
set(gca, 'xticklabel', til, 'FontSize', 10, 'FontWeight', 'bold');
set(gca, 'XTickLabelRotation', 45);
xlabel('k [{{2\pi}/{a}}]', 'Fontsize', 12, 'FontWeight', 'bold');
ylabel('E [eV]', 'Fontsize', 12, 'FontWeight', 'bold');
grid on;
title(semiconductor, 'FontSize', 12, 'FontWeight', 'bold');
waitfor(graph_lines);
toc


