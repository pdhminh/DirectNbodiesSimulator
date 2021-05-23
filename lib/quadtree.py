class Node():
    def __init__(self):
        self.children = [None, None, None, None]

    def insert(self, payload, xo, yo, w, h):
        x = payload.get_x
        y = payload.get_y
        index = self.find_quadnum(x, y, xo, yo, w, h)
        target = self.children[index]
        if target is None:
            target = payload
        elif isinstance(target, Node):
            target.insert(payload, xo, yo, w / 2, h / 2)
        else:
            extra = target
            if index == 0:
                xo = xo - w / 2
                yo = yo + h / 2
            if index == 1:
                xo = xo + w / 2
                yo = yo + h / 2
            if index == 2:
                xo = xo + w / 2
                yo = yo - h / 2
            if index == 3:
                xo = xo - w / 2
                yo = yo - h / 2
            target = Node()
            target.inset(extra, xo, yo, w / 2, h / 2)
            target.insert(payload, xo, yo, w / 2, h / 2)

    def find_quadnum(self, x, y, xo, yo, w, h):
        if xo - w / 2 <= x < xo:
            if y >= yo:
                return 0
            else:
                return 3
        elif xo <= x < xo + w / 2:
            if y >= yo:
                return 1
            else:
                return 2
