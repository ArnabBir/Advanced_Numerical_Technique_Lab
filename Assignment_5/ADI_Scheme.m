k=1;
d = 1.0/4.0;
r = 1.0/2.0;

rx = r;
ry = r;

dx = d;
dy = d;

dt = r*d*d/k;

% Dividing time mesh into two part
dt = dt/2;

% Initial Condition
i_u = 0;
% Boundary Condition
b_u = @(t,x,y)cos(pi*x/2.0)* cos(pi*y/2.0);

u0 = 0;
un = 0;

x0 = 0;
y0 = 0;
t0 = 0;

xn = 1;
yn = 1;
tn = 0.4;

x = x0:dx:xn;
y = y0:dy:yn;
t = t0:dt:tn;

i = length(x) - 1;
j = length(y) - 1;
n = length(t) - 1;

u = zeros(i+1, j+1 , n+1);

% Setting Initial Condition
u(:,:, 1) = 0;

% Setting Boundary Condition
for a = 1:i+1,
    for c = 2:n+1,
        u(a,1,c) = b_u(t(c),x(a),y(1));
        u(a,j+1,c) = b_u(t(c),x(a),y(j+1));
    end
end
for b = 1:j+1,
    for c = 2:n+1,
        u(1,b,c) = b_u(t(c),x(1),y(b));
        u(i+1,b,c) = b_u(t(c),x(i+1),y(b));
    end
end

% Solution

ai_x = @(x)( rx);
bi_x = @(x)( -2*(ry+1) );
ci_x = @(x)( rx );
di_x = @(pmesh, x, y)( -2*pmesh(x,y) - rx*(dx/dy)*(dx/dy)*(pmesh(x,y+1) - 2*pmesh(x,y) + pmesh(x+1,y)) );

ai_y = @(y)( ry );
bi_y = @(y)( -2*(ry+1) );
ci_y = @(y)( ry );
di_y = @(pmesh, x, y)( -r*pmesh(x-1, y) + 2*(r-1)*pmesh(x, y) - r*pmesh(x+1, y) );

disp(n)
for a = 2:n+1,

    if mod(a,2)==0,
        for b = 2:j,
            mat_A = zeros(i-1, i-1);
            vec_b = zeros(i-1, 1);

            % Preparing trigonal matrix
            for g=1:(i),
                mat_A(g,g) = bi_x(NaN);
                mat_A(g,g+1) = ci_x(NaN);
                mat_A(g+1,g) = ai_x(NaN);

                if (g~=i),
                    vec_b(g) = di_x(u(:,:,a-1),g+1, b );
                end
            end

            vec_b(1) = vec_b(1) - rx*u(1,b,a-1); % Check this
            vec_b(i-1) = vec_b(i-1) - rx*u(i,b,a-1);

            mat_A = mat_A(1:i-1, 1:i-1);

            vec_b = vec_b(1:i-1);

            u(2:i, b, a) = thomas_algorithm(mat_A, vec_b);
        end
    
    else
        
        for b = 2:i,
            mat_A = zeros(j-1, j-1);
            vec_b = zeros(j-1, 1);

            % Preparing trigonal matrix
            for g=1:(j),
                mat_A(g,g) = bi_y(NaN);
                mat_A(g,g+1) = ci_y(NaN);
                mat_A(g+1,g) = ai_y(NaN);

                if (g~=j),
                    vec_b(g) = di_y(u(:,:,a-1),b , g+1 );
                end
            end

            vec_b(1) = vec_b(1) - ry*u(b,1,a-1); % Check this
            vec_b(j-1) = vec_b(j-1) - ry*u(b,j,a-1);

            mat_A = mat_A(1:j-1, 1:j-1);

            vec_b = vec_b(1:j-1);

            u(b, 2:j, a) = thomas_algorithm(mat_A, vec_b);
        end 
    end 
end
u

display(x, 'X');
display(y, 'Y');
mesh(x,y,u(:,:,26))
% xlabel('Time')
% ylabel('X')
zlabel('U(x) (Solution)')

% plot(x, y, u(x,y,n));

% plot(X, Y, X, Y, '*', X, Y, 'r');
 xlabel('X(Domain)');
 ylabel('Y(Domain)');
 title('3D - Plot');
