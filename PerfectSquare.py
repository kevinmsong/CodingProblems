class Solution:
    def __init__(self):
        self.product = None
        self.is_perfect_square = False
        self.previous_num = -1
        self.previous_previous_num = -2
    
    def get_perfect_square(self, num):
        print(num, self.product)
        if int(num) ** 2 == self.product:
            self.is_perfect_square = True
            return
        
        if num == self.previous_num or num == self.previous_previous_num:
            return
        
        else:
            self.previous_previous_num = self.previous_num
            self.previous_num = num
            self.get_perfect_square(int(0.5 * (num + self.product / num)))
    
    def isPerfectSquare(self, num: int) -> bool:
        self.product = num
        self.get_perfect_square(1)
        return self.is_perfect_square
