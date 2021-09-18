lambda = function(tau, choice){
  
  if(choice=='sqrt'){
    
    lambda_choice = expression(sqrt(tau))  
    
  } else if(choice=='cubic_root'){
    
    lambda_choice = expression((tau)^(1/3))
    
  } else if(choice=='6th_root'){
    
    lambda_choice = expression((tau)^(1/6))
    
  } else if(choice=='cos'){
    
    lambda_choice = expression(cos(tau))
    
  } else if(choice=='sigmoid'){
    
    lambda_choice = expression(1/(1+exp(-tau)))
    
  } else if(choice=='hacovercos'){
    
    lambda_choice = expression((sin(tau) + 1)/2)
    
  } else if(choice=='sum_of_roots'){
    
    lambda_choice = expression((tau)^(1/3) + (tau)^(1/2))
  }
  
  return(lambda_choice)
  
} 