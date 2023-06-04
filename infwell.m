clc, clearvars,clf,cla
Lx = 20;
Ly = 21;
Lz = 10;
nevs = 10; % number of eigenpairs to find
nvecs = 8; % number of eigenvectors to plot

meff = 0.04; % units of electron mass
a = 0.5;   % nm
t = 0.0380998; % [eV/(nm^2 me)unit prefactor to get eV

H=(t/a^2).*hamilt(Lx, Ly, Lz);
Hp= (t/a^2).*hamiltperiodic(Lx, Ly,Lz);
identity = speye(size(H));

[~,D] = eigs(H,identity,nevs,"smallestabs");
[V,D1] = eigs(Hp,identity,nevs*2,"smallestabs");
D = diag(D);
D1all = D1;
D1 = unique(diag(D1));
disp("fixed, D")
disp(D)
disp("periodic, D1")
disp(D1)
%% plot eigen energies
ncols = 2;
nrows = ceil(nvecs/ncols)+1;
subplot(nrows, ncols,[1,2])
l1 = plot(nan, nan, "b-");
hold on
l2 = plot(nan, nan, "r--");
hold off
legend([l1,l2], ["fixed", "periodic"])
xline(D, "b-","LineWidth",1,'HandleVisibility','off');
xline(D1, "r--", "LineWidth",1,'HandleVisibility','off');
xlabel("Eigen energies (eV)");


%% plot eigenvectors
[Y,X,Z] = meshgrid(1:Ly, 1:Lx, 1:Lz);

for k = [1:nvecs]
    e1 = reshape(V(:,k), Lx, Ly, Lz);
    e1 = smooth3(e1);

    subplot(nrows,ncols,k+ncols)
    slice(Y,X,Z,e1,[Ly/2],[Lx/2],[Lz/2], "cubic")

    colorbar()
    view(3)
    axis vis3d
    xlabel("y")
    ylabel("x")
    zlabel("Z")
    title(D1all(k,k))
    grid off
    shading flat
end


%% functions
function f = index(nx,ny,nz)
    f = @(x, y, z) x + (y-1).*nx + (z-1).*(ny*nx);
end

function H = hamilt(nx, ny, nz)
    [X,Y,Z] = meshgrid(1:nx, 1:ny, 1:nz);
    ind     = index(nx,ny,nz);
    X = reshape(X,  1, []);
    Y = reshape(Y, 1, []);
    Z = reshape(Z, 1, []);
    N = size(X,2);

    H = sparse(N,N);
    for k = 1:N
        x = X(k);
        y = Y(k);
        z = Z(k);
        row = ind(x,y,z);

        if (x ~= 1)
            H(row,ind(x-1,y,z)) = -1;
        end
        if (x ~= nx)
            H(row,ind(x+1,y,z)) = -1;
        end
        if (y ~= 1)
            H(row,ind(x,y-1,z)) = -1;
        end
        if (y ~= ny)
            H(row,ind(x,y+1,z)) = -1;
        end
        if (z ~= 1)
            H(row,ind(x,y,z-1)) = -1;
        end
        if (z ~= nz)
            H(row,ind(x,y,z+1)) = -1;
        end
        H(row,ind(x,y,z)) = 6;
    end
end

function H = hamiltperiodic(nx, ny, nz)
    [X,Y,Z] = meshgrid(1:nx, 1:ny, 1:nz);
    ind     = index(nx,ny,nz);
    X = reshape(X,  1, []);
    Y = reshape(Y, 1, []);
    Z = reshape(Z, 1, []);
    N = size(X,2);

    H = sparse(N,N);
    for k = 1:N
        x = X(k);
        y = Y(k);
        z = Z(k);
        row = ind(x,y,z);

        if (x ~= 1)
            H(row,ind(x-1,y,z)) = -1;
        else
            H(row,ind(nx,y,z)) = -1;
        end

        if (x ~= nx)
            H(row,ind(x+1,y,z)) = -1;
        else
            H(row,ind(1,y,z)) = -1;
        end

        if (y ~= 1)
            H(row,ind(x,y-1,z)) = -1;
        else
            H(row,ind(x,ny,z)) = -1;
        end

        if (y ~= ny)
            H(row,ind(x,y+1,z)) = -1;
        else
            H(row,ind(x,1,z)) = -1;
        end

        if (z ~= 1)
            H(row,ind(x,y,z-1)) = -1;
        else
            H(row,ind(x,y,nz)) = -1;
        end

        if (z ~= nz)
            H(row,ind(x,y,z+1)) = -1;
        else
            H(row,ind(x,y,1)) = -1;
        end
        H(row,ind(x,y,z)) = 6;
    end
end

function D = dlap(n)
    % add -1 values on the left
    n2=n-2;
    range = 1:n2;
    i = range;
    j = range;
    v = -ones(1,n2);
    % add -1 values on the right
    i = [i, range];
    j = [j, (range + 2)];
    v = [v, -ones(1,n2)];
    % add 2 values in the middle
    i = [i, range];
    j = [j, (range + 1)];
    v = [v, 2.*ones(1,n2)];
    
    D = sparse(i+1, j, v, n, n);
end

