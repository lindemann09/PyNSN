import pynsn

r = pynsn.Rectangle((10, 20), (23, 78))
d = pynsn.Dot([10, 20], 78)

print(d.colour.value)
print(d)
d.xy = [100, 100]
print(r.left_top)
print(r.right_bottom)

a = np.array(
    [
        [10, 12.43],
        [10, 12.43],
        [np.nan, np.nan],
        [10, 12.43],
        [np.nan, n.nan],
        [10, 12.43],
        [10, 12.43],
    ]
)
