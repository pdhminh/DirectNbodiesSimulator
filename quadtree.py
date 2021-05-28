# from pprint import pprint
from decimal import Decimal

class Node():
    def __init__(self):
        self.children = [None, None, None, None]
        self.gx = 0
        self.gy = 0
        self.mass = 0

    def get_x(self):
        return self.gx

    def get_y(self):
        return self.gy

    def insert(self, payload, xo, yo, w, h): 
        x = payload.get_x()
        y = payload.get_y()
        print('(', x,',', y,')', '(', xo,',', yo,')', w, h)
        index = self.find_quadnum(x, y, xo, yo, w, h)
        target = self.children[index]
        if target is None:
            self.children[index] = payload
        else:
            if isinstance(target, Node):
                xo, yo, w, h = self.updatexoyo(xo, yo, w, h, index)
                target.insert(payload, xo, yo, w, h)
            else:
                xo, yo, w, h = self.updatexoyo(xo, yo, w, h, index)
                self.children[index] = Node()
                self.children[index].insert(target, xo, yo, w, h)
                self.children[index].insert(payload, xo, yo, w, h)
                # target.cal_approxiamtions(target  )

    def find_quadnum(self, x, y, xo, yo, w, h):
        # x = Decimal(x)
        # y = Decimal(y)
        # xo = Decimal(xo)
        # yo = Decimal(yo)
        # w = Decimal(w)
        # h = Decimal(h)
        if xo - (w / 2) <= Decimal(x) < xo:
            print('left')
            if Decimal(y) >= yo:
                return 0
            else:
                return 3
        if xo <= Decimal(x) <= xo + (w / 2):
            print('right')
            if Decimal(y) >= yo:
                return 1
            else:
                return 2

    def cal_approxiamtions(self, payload):
        if not isinstance(payload, Node):
            return
        for children in payload.children:
            self.cal_approxiamtions(children)
        M = 0
        gx = 0
        gy = 0
        for children in payload.children:
            if children is not None:
                M = M + children.mass
                gx = gx + children.mass * children.get_x()
                gy = gy + children.mass * children.get_y()
        payload.gx = gx / M
        payload.gy = gy / M
        payload.mass = M

    def remove(self, payload):
        for i in range(0, 4, 1):
            if self.children[i] == payload:
                self.children[i] = None
            if isinstance(self.children[i], Node):
                self.children[i].remove(payload)
            if isinstance(self.children[i], Node) and self.children[i].children == [None, None, None, None]:
                self.children[i] = None
        return self

    def updatexoyo(self, xo, yo, w, h, index):
        if index == 0:
            xo = xo - w / 4
            yo = yo + h / 4
        if index == 1:
            xo = xo + w / 4
            yo = yo + h / 4
        if index == 2:
            xo = xo + w / 4
            yo = yo - h / 4
        if index == 3:
            xo = xo - w / 4
            yo = yo - h / 4
        w = w / 2
        h = h / 2 
        return (xo, yo, w, h)
