import pynsn
from tests.shapes_test_picture import shapes_test_picture, Line


d_big = pynsn.Dot((-8, 9), diameter=120, attribute="#00FF00")
da = pynsn.Dot((-20, 70), diameter=60, attribute="#222800")
db = pynsn.Dot((-63, -50), diameter=40, attribute="#000088")

ra = pynsn.Rectangle((50, 75), size=(40, 40), attribute="#FF0000")
rb = pynsn.Rectangle((-50, -45), size=(40, 40), attribute="#ccFF00")
r_big = pynsn.Rectangle((10, -40), size=(150, 60), attribute="#000088")

l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")
shapes_test_picture((d_big, ra, db, l))


print(d_big.intersects(ra))
