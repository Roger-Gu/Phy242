# do integeration on func(x,y,arg1,..,), from x = a to b, y from gfun(x) to hfun(x)
def dblquad(func, a, b, gfun, hfun, args=(), x_step = 0.01, y_step = 0.01):
    x = a
    result = 0
    while (x < b):
        y = gfun(x)
        y_max = hfun(x)
        while (y < y_max):
            result += (func(x,y,*args) * y_step * x_step)
            y += y_step
        x += x_step
    return result