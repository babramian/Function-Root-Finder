import numpy as np

class Equation():

    equation = None
    true_roots = None
    true_root = None


    bisec_ea = []
    NR_ea = []
    secant_ea = []
    FP_ea = []
    mod_ea = []

    bisec_et = []
    NR_et = []
    secant_et = []
    FP_et = []
    mod_et = []

    bisec_num_iter = None
    NR_num_iter = None
    secant_num_iter = None
    FP_num_iter = None
    mod_num_iter = None

    def __init__(self, equation, true_roots):
        self.equation = equation
        self.true_roots = true_roots
    


    def solve(self, num):
        return self.equation(float(num))


    def get_root(self, roots, guess):
        temp = [(x-guess) for x in roots]
        temp2 = None
        if len(roots) > 1:
            temp2 = roots[temp.index(min(temp))]
            return temp2

        else:
            temp2 = roots[0]
            return temp2


    # Bisection Method
    # Input: Starting values A and B, Maximum error value, i is max iterations
    # Returns: Root of the equation between A and B
    def bisection(self, a, b, err=.01, N=100, true_root=[]):
        A = float(a)
        B = float(b)
        count = 1
        max_error = err

        if self.solve(A) * self.solve(B) > 0:
            print("Initial guesses were incorrect.")

        C = (A + B)/2
       
        if true_root == []:
            true_root = self.get_root(self.true_roots, C)

        if self.solve(A) * self.solve(C) < 0:
            B = C
        else:
            A = C

        approx_error = abs(((A + B)/2 - C)/((A + B)/2))
        true_error = abs((true_root - C)/true_root)
        self.bisec_ea.append(float(approx_error))
        self.bisec_et.append(float(true_error))

        while approx_error > max_error and count <= N:
            if count == N:
                print(f"Maximum iterations reached ({count})")
                break

            C = (A + B)/2
            #print(f'C: {C}\t\tA: {A}\t\tB: {B}')

            if self.solve(A) * self.solve(C) < 0:
                B = C
            else:
                A = C
            
            approx_error = abs(((A + B)/2 - C)/((A + B)/2))
            true_error = abs((true_root - C)/true_root)
            self.bisec_ea.append(float(approx_error))
            self.bisec_et.append(float(true_error))
            count += 1
            
        self.bisec_num_iter = count
        return C
    
    

    def newton_raphson(self, derivative, x0, err=.01, N=100, true_root=[]):
        step = 0
        flag = 1
        ea = err + 1
        x1 = None

        if true_root == []:
            true_root = self.get_root(self.true_roots, x0 - self.solve(x0)/derivative(x0))
        
        while step < N and ea > err:
            if derivative(x0) == 0.0:
                print("Can't divide by zero.")
                break

            x1 = x0 - self.solve(x0)/derivative(x0)
            #print(f'Iteration {step+1},\tx0 = {round(x0, 3)}\tx1 = {round(x1, 3)}\tf(x1) = {round(self.solve(x1), 3)}')
            
            ea = abs((x1 - x0)/x1)
            et = abs((true_root - x1)/true_root)
            self.NR_ea.append(float(ea))
            self.NR_et.append(float(et))

            x0 = x1
            step += 1

        self.NR_num_iter = step

        if step <= N:
            print()
            print(f"(NR) Root is: {x1}")
            print()
            return x1

        else:
            print("Not Convergent")



    def secant(self, x0, x1, err=0.01, N=100, true_root=[]):
        step = 1
        ea = err + 1
        x2 = None

        if true_root == []:
            true_root = self.get_root(self.true_roots, x1 - (x1 - x0) * self.solve(x1)/(self.solve(x1) - self.solve(x0)))

        while step < N and ea > err:
            if self.solve(x0) == self.solve(x1):
                print("Mathematical Error")
                break
            
            print(f'Iteration {step},\tx0 = {round(x0, 3)}\tx1 = {round(x1, 3)}\tf(x1) = {round(self.solve(x1), 3)}')

            ea = abs((x1 - x0)/x1)
            et = abs((true_root - x1)/true_root)
            self.secant_ea.append(float(ea))
            self.secant_et.append(float(et))

            x2 = x1 - (x1 - x0) * self.solve(x1)/(self.solve(x1) - self.solve(x0))
            x0 = x1
            x1 = x2

            step += 1

        print(f'Iteration {step},\tx0 = {round(x0, 3)}\tx1 = {round(x1, 3)}\tf(x1) = {round(self.solve(x1), 3)}')
        ea = abs((x1 - x0)/x1)
        et = abs((true_root - x1)/true_root)
        self.secant_ea.append(float(ea))
        self.secant_et.append(float(et))
        self.secant_num_iter = step

        if step <= N:
            print()
            print(f"(Secant) Root is: {x1}")
            print()
            return x1

        else:
            print("Not Convergent")
        
    

    def false_position(self, a, b, err=.01, N=100, true_root=[]):
        step = 0
        ea = err + 1
        #c = a - ((a - b) * self.solve(a))/(self.solve(a) - self.solve(b))

        #if self.solve(a) * self.solve(c) < 0:
        #    b = c
        #else:
        #    a = c

        if self.solve(a) * self.solve(b) > 0:
            print("Incorrect initial guesses")
            print(f"Pick a new value for A: ")
            a = float(input())
            print(f"Pick a new value for B: ")
            b = float(input())

        if true_root == []:
            true_root = self.get_root(self.true_roots, a - ((a - b) * self.solve(a))/(self.solve(a) - self.solve(b)))

        while ea > err and step < N:
            print(f'Iteration {step+1}')

            c = a - ((a - b) * self.solve(a))/(self.solve(a) - self.solve(b))
            print(f'C: {c}')

            if self.solve(a) * self.solve(c) < 0:
                b = c
            else:
                a = c

            ea = abs(((a - ((a - b) * self.solve(a))/(self.solve(a) - self.solve(b))) - c)
                    /(a - ((a - b) * self.solve(a))/(self.solve(a) - self.solve(b))))
            et = abs((true_root - (c))/true_root)
            self.FP_ea.append(float(ea))
            self.FP_et.append(float(et))

            step += 1

        self.FP_num_iter = step
        if step <= N:
            print()
            print(f"(FP) Root is: {c}")
            print()
            return c

        else:
            print("Not Convergent")
            

    
    def reset_values(self):
        self.bisec_ea = []
        self.NR_ea = []
        self.secant_ea = []
        self.FP_ea = []
        self.mod_ea = []

        self.bisec_et = []
        self.NR_et = []
        self.secant_et = []
        self.FP_et = []
        self.mod_et = []

        self.bisec_num_iter = None
        self.NR_num_iter = None
        self.secant_num_iter = None
        self.FP_num_iter = None
        self.mod_num_iter = None



    def get_errors(self):
        return (np.array(self.bisec_ea, dtype=float), np.array(self.NR_ea, dtype=float), 
                np.array(self.secant_ea, dtype=float), np.array(self.FP_ea, dtype=float),
                np.array(self.mod_ea, dtype=float),
                np.array(self.bisec_et, dtype=float), np.array(self.NR_et, dtype=float), 
                np.array(self.secant_et, dtype=float), np.array(self.FP_et, dtype=float),
                np.array(self.mod_et, dtype=float))
        
    



    def get_max_iterations(self):
        return np.array([self.bisec_num_iter, 
                            self.NR_num_iter,
                            self.secant_num_iter,
                            self.FP_num_iter,
                            self.mod_num_iter], dtype=float)

    


if __name__ == "__main__":
    # True Roots of fa = .365, 1.922, 3.563 (According to Desmos)
    fa = Equation((lambda x : 2 * (x**3) - 11.7 * (x**2) + 17.7 * x - 5), [.365, 1.922, 2.563])
    fa_derivative = lambda x : 6 * (x**2) - 23.4 * x + 17.7

    # True Root of fb = 126.632 (According to Desmos)
    fb = Equation((lambda x : x + 10 - x * np.cosh(50/x)), [126.632])
    fb_derivative = lambda x : 1 + (50 * np.sinh(50/x))/x - np.cosh(50/x)


    bis_a1 = fa.bisection(0, 1, true_root=.365)
    NR_a1 = fa.newton_raphson(fa_derivative, -5, true_root=.365)
    secant_a1 = fa.secant(-10, -5, true_root=.365)
    FP_a1 = fa.false_position(-2, 1, true_root=.365)
    bis_ea_a1, NR_ea_a1, secant_ea_a1, FP_ea_a1, mod_ea_a1, bis_et_a1, NR_et_a1, secant_et_a1, FP_et_a1, mod_er_a1 = fa.get_errors()
    bis_count_a1, NR_count_a1, secant_count_a1, FP_count_a1, mod_count_a1 = fa.get_max_iterations()
    fa.reset_values()

    bis_a2 = fa.bisection(1, 3, true_root=1.922)
    NR_a2 = fa.newton_raphson(fa_derivative, 2.5, true_root=1.922)
    secant_a2 = fa.secant(.5, 1.5, true_root=1.922)
    FP_a2 = fa.false_position(1.4, 2.3, true_root=1.922)
    bis_ea_a2, NR_ea_a2, secant_ea_a2, FP_ea_a2, mod_ea_a2, bis_et_a2, NR_et_a2, secant_et_a2, FP_et_a2, mod_et_a2 = fa.get_errors()
    bis_count_a2, NR_count_a2, secant_count_a2, FP_count_a2, mod_count_a2 = fa.get_max_iterations()
    fa.reset_values()
    
    bis_a3 = fa.bisection(3, 4, true_root=3.568)
    NR_a3 = fa.newton_raphson(fa_derivative, 10, true_root=3.568)
    secant_a3 = fa.secant(40, 25, true_root=3.568)
    FP_a3 = fa.false_position(3, 5, true_root=3.568)
    bis_ea_a3, NR_ea_a3, secant_ea_a3, FP_ea_a3, mod_ea_a3, bis_et_a3, NR_et_a3, secant_et_a3, FP_et_a3, mod_et_a3 = fa.get_errors()
    bis_count_a3, NR_count_a3, secant_count_a3, FP_count_a3, mod_count_a3 = fa.get_max_iterations()
    fa.reset_values()
    fb.reset_values()

    bis_b = fb.bisection(123, 130, true_root=126.632)
    secant_b = fb.secant(100, 110, true_root=126.632)
    NR_b = fb.newton_raphson(fb_derivative, 120, true_root=126.632)
    FP_b = fb.false_position(120, 250, true_root=126.632)
    bis_ea_b, NR_ea_b, secant_ea_b, FP_ea_b, mod_ea_b, bis_et_b, NR_et_b, secant_et_b, FP_et_b, mod_et_b = fb.get_errors()
    bis_count_b, NR_count_b, secant_count_b, FP_count_b, mod_count_b = fb.get_max_iterations()

    print()
    print(f'Bisections of function (a) {round(bis_a1, 3)}, {round(bis_a2, 3)}, {round(bis_a3, 3)}')
    print(f'f(bisections): {round(fa.solve(bis_a1), 3)}, {round(fa.solve(bis_a2), 3)}, {round(fa.solve(bis_a3), 3)}')
    print(f'Bisections of function (a) {round(bis_b, 3)}')
    print(f'f({bis_b}): {round(fb.solve(bis_b), 3)}')
    print()
    print(f'A1: {FP_ea_a1, FP_et_a1}')
    print(f'A2: {FP_ea_a2, FP_et_a2}')
    print(f'A3: {FP_ea_a3, FP_et_a3}')
    print(f'B: {FP_ea_b, FP_et_b}')
    print(f'FP counts: {FP_count_a1, FP_count_a2, FP_count_a3, FP_count_b}')
    
    


input()