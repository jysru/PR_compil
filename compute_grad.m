function grad = compute_grad(z, u, y, Amatrix,tt)
    ayz = Amatrix(tt,:)*z;
    yz = sqrt(abs(ayz).^2+u^2);
    grad = Amatrix(tt,:)'*(((yz-y(tt))./yz).*ayz);
end
      
