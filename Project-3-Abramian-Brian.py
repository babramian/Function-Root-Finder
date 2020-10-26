import numpy as np

class Equation():

    equation = None

    def __init__(self, equation):
        self.equation = equation
    
    def solve(self, num):
        return self.equation(float(num))

    # Bisection Method
    # Input: Starting values A and B, Maximum error value
    # Returns: Root of the equation between A and B
    def bisection(self, a, b, err):
        A = float(a)
        B = float(b)
        error = err

        if (self.solve(A) * self.solve(B) > 0):
            print("Initial guesses were incorrect.")

        C = (A + B)/2

        if (self.solve(A) * self.solve(C) < 0):
            B = C
        else:
            A = C

        while (abs(((A + B)/2 - C)/((A + B)/2)) > error):
            C = (A + B)/2
            #print(f'C: {C}\t\tA: {A}\t\tB: {B}')

            if (self.solve(A) * self.solve(C) < 0):
                B = C
            else:
                A = C
        
        return C
        


if __name__ == "__main__":
    # True Roots of fa = .365, 1.922, 3.563 (According to Desmos)
    fa = Equation((lambda x : 2 * (x**3) - 11.7 * (x**2) + 17.7 * x - 5))
    # True Root of fb = 126.632 (According to Desmos)
    fb = Equation((lambda x : x + 10 - x * np.cosh(50/x)))

    print(fa.solve(0))
    print(fa.solve(1))
    print(fa.solve(.365))
    bis_a1 = fa.bisection(0, 1, .001)
    bis_a2 = fa.bisection(1, 3, .001)
    bis_a3 = fa.bisection(3, 4, .001)
    print()
    print(f'Bisections of function (a) {round(bis_a1, 3)}, {round(bis_a2, 3)}, {round(bis_a3, 3)}')
    print(f'f(bisections): {round(fa.solve(bis_a1), 3)}, {round(fa.solve(bis_a2), 3)}, {round(fa.solve(bis_a3), 3)}')

    
input()