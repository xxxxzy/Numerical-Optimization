import numpy as np
import scipy
import scipy.linalg
import matplotlib.pyplot as plt


x =[1,1]
x1 = []
y = []
deltap = []
iteration = 0
r = 10**(-8)

while 1 > 0:
    object_fx = lambda x : np.log1p((x[0]**2 + 1000*x[1]**2)/e)  #to input the objective function
   #to input the initial seaching point
    e = 10**(-16)

    epsilon =10**(-4) 

    delta = 1



    x_k = x

    fx = object_fx

    def G(x_k):
        
        alpha1 = 1000
        epsilon = 10**(-16)
    
        gra = np.array([2*x_k[0]/(epsilon+x_k[0]**2+alpha1*x_k[1]**2),2*alpha1*x_k[1]/(epsilon+x_k[0]**2+alpha1*x_k[1]**2)])
#        print(gra)

        return gra


    hess = np.array([[5,-1],[-4,4]])
    


    if np.all(np.linalg.eigvals(hess)) > 0:
    
        p_0 = np.dot(-np.linalg.inv(hess),G(x_k))
        p = p_0
    
        if np.linalg.norm(p_0) <= delta**2:
        
        
            lamda_l = 1
            iteration = 0
        
            while 1>0:
            
                RR = hess + lamda_l * np.diag([1,1])
                RT = scipy.linalg.cholesky(RR, lower=True)
            
            
                p_l = -np.dot(G(x_k),np.linalg.inv(RR))
                q_l = np.dot(p_l,np.linalg.inv(RT))

            
                p = np.linalg.norm(p_l)
                q = np.linalg.norm(q_l)
              
                p_lamda = -np.dot(G(x_k),np.linalg.inv(hess + lamda_l * np.diag([1,1])))
          
            
                diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*(np.linalg.norm(p_l)-delta)/delta
            
                iteration = iteration + 1
                if iteration == 2 or not np.all(hess + (lamda_l + diff)* np.diag([1,1])) > 0:
                    lamda_0 = 1
                    lamda_l = lamda_l
#                print(lamda_l)
                    break
                else:
                    lamda_l = lamda_l + diff
                
        else:
        
            lamda_l = 0
            while 1>0:
            
                RR = hess +lamda_l*np.diag([1,1])
                RT = scipy.linalg.cholesky(RR, lower=True)
            
                p_lamda = -np.dot(G(x_k),np.linalg.inv(hess + lamda_l * np.diag([1,1])))
            
                p_l = -np.dot(G(x_k),np.linalg.inv(RR))
                q_l = np.dot(p_l,np.linalg.inv(RT))
              
                if np.linalg.norm(p_lamda) <= delta:
                
                    lamda_0 = 0
                    lamda_l = lamda_l
                
                    break
                else:
                    diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*2*(np.linalg.norm(p_l)-delta)/delta

                    lamda_l = lamda_l + diff
                
#            print(lamda_l)
                
                
                


    else:
    
        lamda_l = 0
        while 1>0:
            RR =hess+lamda_l*np.diag([1,1])
#            print(RR)
            RT = scipy.linalg.cholesky(RR, lower=True)
        
            p_l = -np.dot(G(x_k),np.linalg.inv(RR))
            q_l = np.dot(p_l,np.linalg.inv(RT))
        
            equa = hess + lamda_l * np.diag([1,1])
            p_lamda = -np.dot(G(x_k),np.linalg.inv(equa))
            
            if np.linalg.norm(p_lamda) > delta and np.all(np.linalg.eigvals(equa)) > 0:
            
                c = lamda_l
                lamda_0 = lamda_l
            
#            print("lamda_0",lamda_0)
               
                break
            else:
            
                diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*2*(np.linalg.norm(p_l)-delta)/delta
            
                lamda_l = lamda_l + diff
            
            
            
        while 1>0:
            
        
            RR = hess+lamda_l*np.diag([1,1])
#            print(RR)
            RT = scipy.linalg.cholesky(RR, lower=True)
        
            p_l = -np.dot(G(x_k),np.linalg.inv(RR))
            q_l = np.dot(p_l,np.linalg.inv(RT))
        
            equa = hess + lamda_l * np.diag([1,1])
            p_lamda = -np.dot(G(x_k),np.linalg.inv(equa))
            
            if np.linalg.norm(p_lamda) < delta and lamda_l > c:
            
#            print("lamda_l",lamda_l)
               
                break
            else:
            
                diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*2*(np.linalg.norm(p_l)-delta)/delta
            
                lamda_l = lamda_l + diff

    while 1>0:

        if lamda_l - lamda_0 <=10**(-2):
        
            lamda = (lamda_l + lamda_0)/2

    
        else:
        
            lamda = (lamda_0 + lamda_l)/2
            equa = hess + lamda * np.diag([1,1])
            p_lamda = -np.dot(G(x_k),np.linalg.inv(equa))
        
            if np.linalg.norm(p-lamda) > delta:
            
                lamda_0 = lamda
            else:
                lamda_1 = lamda
                
        iteration =iteration + 1
        x1.append(iteration)
        
        norm = np.log(np.linalg.norm(G(x_k)))
        y.append(norm)
        
        deltap.append(delta)
        
                      
        if np.dot(G(x_k).T,G(x_k)) <= epsilon or iteration == 500:
            
            break
        else:
            
            
            p_lamda = -np.dot(G(x_k),np.linalg.inv(hess + lamda * np.diag([1,1])))
            mx = fx(x_k) + np.dot(G(x_k),p_lamda.T) +0.5*np.dot(np.dot(p_lamda.T,hess),p_lamda)
        

            
            rho_l = (fx(x_k)-fx(x_k+p_lamda))/(fx(x_k)-mx)
        

    
            if rho_l < 0.25:
                delta = 0.25*delta
             
            else:
            
                if rho_l >0.75 and np.linalg.norm(p_lamda) == delta:
                    delta = min(2*delta,10)
                
                else:
                    delta = delta
                    
            
            if rho_l > 0:
            
                x_k = x_k + p_lamda
                y_k = G(x_k) - G(x_k-p_lamda)
#                print(x_k)
        
            else:
                x_k = x_k
                y_k = 0
#            print(x_k)
#            print(np.dot(G(x_k).T,G(x_k)))
            
            sec = np.mat(y_k).T-np.dot(hess,np.mat(p_lamda).T)    
            
            if abs(np.dot(np.mat(p_lamda),sec)) >= r*np.linalg.norm(p_lamda)*np.linalg.norm(sec):
                
                
                hess = hess + np.dot(sec,sec.T)/np.dot(sec.T,p_lamda)
                hess = np.array(hess)
#                print('hess',hess)
#                print(type(hess))
                
            else:
                    
                hess = hess
                
#            print(hess)
    print(x_k)            
    print(iteration)
    plt.plot(x1,y)
    plt.legend(labels = ['f3'],loc = 'upper right')
    plt.xlabel("Number of iterations")
    plt.ylabel("log||g||")  
    plt.show()
                

            
    break
            