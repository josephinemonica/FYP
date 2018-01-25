import scipy.optimize as optimize

def f(params):
    # print(params)  # <-- you'll see that params is a NumPy array
    #a, b, c = params # <-- for readability you may wish to assign names to the component variables
    #return a**2 + b**2 + c**2
    a= params # <-- for readability you may wish to assign names to the component variables
    return a**2 
initial_guess = [1]
#b=((-1.5, 1.5), (-1.5, 1.5), (-1.5, 1.5))
b=[(-1.5, 1.5)]
print(len(initial_guess))
print(len(b))
result = optimize.minimize(f, initial_guess,method="L-BFGS-B",bounds=b)
if result.success:
    fitted_params = result.x
    print(fitted_params)
else:
    raise ValueError(result.message)
