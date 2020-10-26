import numpy as np

class Equation():

    equation = None
    true_roots = None
    true_root = None


    bisec_ea = []
    NR_ea = []

    bisec_et = []
    NR_et = []

    bisec_num_iter = None
    NR_num_iter = None

    def __init__(self, equation, true_roots):
        self.equation = equation
        self.true_roots = true_roots
    
    def solve(self, num):
        return self.equation(float(num))

    # Bisection Method
    # Input: Starting values A and B, Maximum error value, i is max iterations
    # Returns: Root of the equation between A and B
    def bisection(self, a, b, err=.01, i=100):
        A = float(a)
        B = float(b)
        count = 1
        max_error = err

        if self.solve(A) * self.solve(B) > 0:
            print("Initial guesses were incorrect.")

        C = (A + B)/2
       
        temp = [(x-C) for x in self.true_roots]
        if len(self.true_roots) > 1:
            self.true_root = self.true_roots[temp.index(min(temp))]

        else:
            self.true_root = self.true_roots[0]

        if self.solve(A) * self.solve(C) < 0:
            B = C
        else:
            A = C

        approx_error = abs(((A + B)/2 - C)/((A + B)/2))
        true_error = abs((self.true_root - C)/self.true_root)
        self.bisec_ea.append(float(approx_error))
        self.bisec_et.append(float(true_error))

        while approx_error > max_error and count <= i:
            if count == i:
                print(f"Maximum iterations reached ({count})")
                break

            C = (A + B)/2
            #print(f'C: {C}\t\tA: {A}\t\tB: {B}')

            if self.solve(A) * self.solve(C) < 0:
                B = C
            else:
                A = C
            
            approx_error = abs(((A + B)/2 - C)/((A + B)/2))
            true_error = abs((self.true_root - C)/self.true_root)
            self.bisec_ea.append(float(approx_error))
            self.bisec_et.append(float(true_error))
            count += 1
            
        self.bisec_num_iter = count
        return C
    
    

    def newton_raphson(self, derivative, x0, err=.01, i=100):
        step = 0
        flag = 1
        ea = err + 1
        x1 = None
        
        while step < i and ea > err:
            if derivative(x0) == 0.0:
                print("Can't divide by zero.")
                break

            x1 = x0 - self.solve(x0)/derivative(x0)
            print(f'Iteration {step+1},\tx0 = {round(x0, 3)}\tx1 = {round(x1, 3)}\tf(x1) = {round(self.solve(x1), 3)}')
            
            ea = abs((x1 - x0)/x1)
            et = abs((self.true_root - x1)/self.true_root)
            self.NR_ea.append(float(ea))
            self.NR_et.append(float(et))

            x0 = x1
            step += 1

        self.NR_num_iter = step

        if step <= i:
            print(f"Root is: {x1}")
            return x1

        else:
            print("Not Convergent")

    
    def reset_values(self):
        self.bisec_ea = []
        self.NR_ea = []

        self.bisec_et = []
        self.NR_et = []

        self.bisec_num_iter = None
        self.NR_num_iter = None


    def get_errors(self):
        return np.array(self.bisec_ea, dtype=float), np.array(self.NR_ea, dtype=float), np.array(self.bisec_et, dtype=float), np.array(self.NR_et, dtype=float)
        
    
    def get_max_iterations(self):
        return np.array([self.bisec_num_iter, self.NR_num_iter], dtype=float)

    


if __name__ == "__main__":
    # True Roots of fa = .365, 1.922, 3.563 (According to Desmos)
    fa = Equation((lambda x : 2 * (x**3) - 11.7 * (x**2) + 17.7 * x - 5), [.365, 1.922, 2.563])
    fa_derivative = lambda x : 6 * (x**2) - 23.4 * x + 17.7

    # True Root of fb = 126.632 (According to Desmos)
    fb = Equation((lambda x : x + 10 - x * np.cosh(50/x)), [126.632])
    fb_derivative = lambda x : 1 + (50 * np.sinh(50/x))/x - np.cosh(50/x)


    bis_a1 = fa.bisection(0, 1)
    NR_a1 = fa.newton_raphson(fa_derivative, -5)
    bis_ea_a1, NR_ea_a1, bid_et_a1, NR_et_a1 = fa.get_errors()
    bis_count_a1, NR_count_a1 = fa.get_max_iterations()
    fa.reset_values()

    bis_a2 = fa.bisection(1, 3)
    NR_a2 = fa.newton_raphson(fa_derivative, 2.5)
    bis_ea_a2, NR_ea_a2, bid_et_a2, NR_et_a2 = fa.get_errors()
    bis_count_a2, NR_count_a2 = fa.get_max_iterations()
    fa.reset_values()
    
    bis_a3 = fa.bisection(3, 4)
    NR_a3 = fa.newton_raphson(fa_derivative, 10)
    bis_ea_a3, NR_ea_a3, bid_et_a3, NR_et_a3 = fa.get_errors()
    bis_count_a3, NR_count_a3 = fa.get_max_iterations()
    fa.reset_values()
    
    fb.reset_values()
    bis_b = fb.bisection(120, 130)
    NR_b = fb.newton_raphson(fb_derivative, 120)
    bis_ea_b, NR_ea_b, bid_et_b, NR_et_b = fb.get_errors()
    bis_count_b, NR_count_b = fb.get_max_iterations()

    print()
    #print(f'Bisections of function (a) {round(bis_a1, 3)}, {round(bis_a2, 3)}, {round(bis_a3, 3)}')
    #print(f'f(bisections): {round(fa.solve(bis_a1), 3)}, {round(fa.solve(bis_a2), 3)}, {round(fa.solve(bis_a3), 3)}')
    #print(f'Bisections of function (a) {round(bis_b, 3)}')
    #print(f'f({bis_b}): {round(fb.solve(bis_b), 3)}')
    #print()
    print(f'A1: {bis_ea_a1, NR_ea_a1, bid_et_a1, NR_et_a1}')
    #print(f'A2: {bis_ea_a2, NR_ea_a2, bid_et_a2, NR_et_a2}')
    #print(f'A3: {bis_ea_a3, NR_ea_a3, bid_et_a3, NR_et_a3}')
    #print(f'B: {bis_ea_b, NR_ea_b, bid_et_b, NR_et_b}')
    print(f'B counts: {bis_count_a1, NR_count_a1}')
    
    


input()